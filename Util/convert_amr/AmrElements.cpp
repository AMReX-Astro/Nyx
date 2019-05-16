//BL_COPYRIGHT_NOTICE

#include <new>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <list>
#include <map>
#include <algorithm>

using std::cout;
using std::cerr;
using std::endl;

#ifndef WIN32
#include <unistd.h>
#endif

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_DataServices.H"
#include "AMReX_Utility.H"
#include "AMReX_FArrayBox.H"

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
    cout << '\n';
    cout << "Usage:" << '\n';
    cout << "    infile  = plotfilename" << '\n';
    cout << "   [-verbose]" << '\n';
    cout << '\n';
    exit(1);

}

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    {

    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }

    std::string infile; 
    pp.get("infile",infile);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Reading " << infile << "...";

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk())  
        DataServices::Dispatch(DataServices::ExitRequest, NULL);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Done reading plot file" << std::endl;

    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel); finestLevel=std::min(finestLevel,amrData.FinestLevel());
    int Nlev = finestLevel + 1;

    if (ParallelDescriptor::IOProcessor())
        std::cout << "... finest level " << finestLevel << std::endl;

    Vector<int> comps;
    if (int nc = pp.countval("comps"))
    {
        comps.resize(nc);
        pp.getarr("comps",comps,0,nc);
    }
    else
    {
        int sComp = 0;
        pp.query("sComp",sComp);
        int nComp = 1;
        pp.query("nComp",nComp);
        BL_ASSERT(sComp+nComp <= amrData.NComp());
        comps.resize(nComp);
        std::cout << "NCOMP NOW " << nComp << std::endl;
        for (int i=0; i<nComp; ++i)
            comps[i] = sComp + i;
    }

    int nComp = comps.size();
    const Vector<string>& plotVarNames=amrData.PlotVarNames();
    Vector<string> inVarNames(nComp);
    Vector<int> destFillComps(nComp);
    for (int i=0; i<nComp; ++i)
    {
        inVarNames[i] = plotVarNames[comps[i]];
        std::cout << "plotVarNames " << plotVarNames[comps[i]] << std::endl;;
        destFillComps[i] = i;
    }

    const int nGrow = 0;

    if (ParallelDescriptor::IOProcessor() && verbose>0)
    {
        cerr << "Will read the following states: ";
        for (int i=0; i<nComp; ++i)
            cerr << " " << amrData.StateNumber(inVarNames[i]) << " (" << inVarNames[i] << ")" ;
        cerr << '\n';
    }

    const Box& probDomain = amrData.ProbDomain()[finestLevel];
    int dir=BL_SPACEDIM-1; pp.query("dir",dir);
    const IntVect lo=probDomain.smallEnd();
    IntVect hi=lo; hi[dir] = probDomain.bigEnd()[dir];
    const Box resBox(lo,hi);
    FArrayBox resFab(resBox,nComp); resFab.setVal(0.);

    int accumFac = 1;
    for (int lev=finestLevel; lev>=0; --lev)
    {
        const BoxArray& ba = amrData.boxArray(lev);
	const DistributionMapping& dm = amrData.DistributionMap(lev);
        MultiFab mf(ba,dm,nComp,nGrow);
        if (ParallelDescriptor::IOProcessor() && verbose>0)
            cerr << "...filling data at level " << lev << endl;
        amrData.FillVar(mf,lev,inVarNames,destFillComps);

	Real bogus_flag=-1e200;
        // Zero covered regions
        if (lev < finestLevel)
        {
            const BoxArray baf = BoxArray(amrData.boxArray(lev+1)).coarsen(amrData.RefRatio()[lev]);

	    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
		FArrayBox& myFab = mf[mfi];
                std::vector< std::pair<int,Box> > isects = baf.intersections(ba[mfi.index()]);
                
		for (int ii = 0; ii < isects.size(); ii++)
		    myFab.setVal(bogus_flag,isects[ii].second,0,nComp);
            }
        }

	///////// Write elements to stdout, to file as ascii, to file as binary
	std::ofstream ofs("fab_ascii.txt", std::ofstream::out);
	std::ofstream ofs_bin("fab_binary", std::ofstream::out);
	std::cout<<"idx\ti\tj\tk\tf(i,j,k,n)"<<std::endl;
	long int idx=0;
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const FArrayBox& fab = mf[mfi];
	    Array4<Real> const& fab_array = mf.array(mfi);
	    int ncomp = mf.nComp();
            const Box& box = mfi.validbox();
	    FArrayBox fab2;
	    fab2.resize(box,ncomp);

            IntVect ivlo = box.smallEnd();
            IntVect ivhi = box.bigEnd(); ivhi[dir] = ivlo[dir];

	    const Dim3 lo = amrex::lbound(box);
	    const Dim3 hi = amrex::ubound(box);
	    const auto len = amrex::length(box);  // length of box

	    for (int n = 0; n < ncomp; ++n) {
	      for (int z = lo.z; z <= hi.z; ++z) {
		for (int y = lo.y; y <= hi.y; ++y) {
		  // Don't use pragma simd
		  //		          AMREX_PRAGMA_SIMD
			    for (int x = lo.x; x <= hi.x; ++x) {
			      //idx assumes ncomp=1, is local to this fab
			      //			      int idx = x+y*len.x+z*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
			      if(fab_array(x,y,z,n) != bogus_flag)
				{
				++idx;
				if(verbose>0)
				  std::cout<<idx<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<fab_array(x,y,z,n) <<std::endl;
				else
				  std::cout<<fab_array(x,y,z,n) <<std::endl;
				}
			    }
		}
	      }
	    }

	    ofs<<fab<<std::endl;
	    fab.writeOn(ofs_bin);

        }

	mf.setVal(0.);

	////////// Read elements from binary file, print to stdout those which don't trip flag
	std::ifstream ifs_bin("fab_binary", std::ofstream::in);
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab2 = mf[mfi];
	    Array4<Real> const& fab_array = mf.array(mfi);
	    int ncomp = mf.nComp();
            const Box& box = mfi.validbox();

            IntVect ivlo = box.smallEnd();
            IntVect ivhi = box.bigEnd(); ivhi[dir] = ivlo[dir];

	    const Dim3 lo = amrex::lbound(box);
	    const Dim3 hi = amrex::ubound(box);

	    fab2.readFrom(ifs_bin);

	    for (int n = 0; n < ncomp; ++n) {
	      for (int z = lo.z; z <= hi.z; ++z) {
		for (int y = lo.y; y <= hi.y; ++y) {
		  // Don't use pragma simd
		  //		          AMREX_PRAGMA_SIMD
			    for (int x = lo.x; x <= hi.x; ++x) {
			      if(fab_array(x,y,z,n) != bogus_flag)
				std::cout<<fab_array(x,y,z,n) <<std::endl;
			    }
		}
	      }
	    }

        }

	//        if (lev>0)

    }

    }
    amrex::Finalize();
    return 0;
}

