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
    cout << "    outfile = outputfilename" << '\n';
    cout << "    comps   = varnumbers" << '\n';
    cout << "   [-verbose]" << '\n';
    cout << '\n';
    exit(1);

}

static Vector<int> findComponents(const AmrData& amrd, const Vector<std::string>& names)
{
    const int n = names.size();
    Vector<int> comps(n);
    for (int i = 0; i < n; i++)
        comps[i] = amrd.StateNumber(names[i]);
    return comps;
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

    bool format_ascii = false;
    //    pp.query("format_ascii",format_ascii);
    if (pp.contains("format_ascii"))
    {
      format_ascii = true;
    }
    //    format_ascii = true;

    std::string infile; 
    std::cout << "Reading " << infile << "...";
    pp.get("infile",infile);

    std::string outfile; 
    pp.get("outfile",outfile);

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
    
    int myRank=0;
#ifdef AMREX_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
#endif

    std::cout<<"finest level: "<<finestLevel<<std::endl;
    std::string filename=outfile+"_"+std::to_string(myRank)+".bin";
    std::string filename_ascii=outfile+"_ascii_"+std::to_string(myRank)+".txt";
    std::ofstream ofs_mpi(filename,std::ios::binary);
    std::ofstream ofs;
    if(format_ascii)
      {
	ofs.open(filename_ascii, std::ofstream::out);
	ofs<<"idx\tx\ty\tz\tdata"<<std::endl;
      }

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

	    const auto dx = amrData.CellSize(lev);
	    RealVect boxSizeH(amrData.ProbHi());
	    if(verbose)
	      {
	    std::cout<<"level"<<lev<<"dx: "<<dx[0]<<"dy: "<<dx[1]<<dx[2]<<std::endl;
	    std::cout<<"Probhi"<<boxSizeH<<std::endl;
	      }
	    /*	    RealBox gridloc = RealBox(ba[mfi.index()],
				      amrData.CellSize(lev),
				      amrData.ProbLo()[lev]);*/

	    /*	    RealBox gridloc = RealBox(ba[mfi.index()],
                                      amrData.CellSize(lev),
                                      amrData.ProbLo());*/

	    const auto len = amrex::length(box);  // length of box

	    if(format_ascii)
	      {
		for (int z = lo.z; z <= hi.z; ++z) {
		  for (int y = lo.y; y <= hi.y; ++y) {
		    // Don't use pragma simd
		    //		          AMREX_PRAGMA_SIMD
		    for (int x = lo.x; x <= hi.x; ++x) {
		      //idx assumes ncomp=1, is local to this fab
		      //			      int idx = x+y*len.x+z*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
		      if(fab_array(x,y,z,0) != bogus_flag)
			{
			  ++idx;
			  if(verbose>0)
			    {
			      ofs<<idx<<"\t"<<x*dx[0]<<"\t"<<y*dx[1]<<"\t"<<z*dx[2]<<"\t"<<dx[0]*dx[1]*dx[2];
			      for (int n = 0; n < comps.size(); ++n) {
				ofs<<"\t"<<fab_array(x,y,z,comps[n]); 
			      }
			      ofs<<std::endl;
			    }
			  else
			    {
			      for (int n = 0; n < comps.size(); ++n) {
				ofs<<"\t"<<fab_array(x,y,z,comps[n]); 
				  }
			      ofs<<std::endl;
			    }
			}
		    }
		  }
		}
	      }

	    for (int z = lo.z; z <= hi.z; ++z) {
	      for (int y = lo.y; y <= hi.y; ++y) {
		// Don't use pragma simd
		//		          AMREX_PRAGMA_SIMD
		for (int x = lo.x; x <= hi.x; ++x) {
		  //idx assumes ncomp=1, is local to this fab
		  //			      int idx = x+y*len.x+z*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
			      if(fab_array(x,y,z,0) != bogus_flag)
				{

				    Real tmpx=x*dx[0];
				    Real tmpy=y*dx[1];
                                    Real tmpz=z*dx[2];
                                    Real tmpvol=dx[2]*dx[0]*dx[1];
				    
				    ofs_mpi.write((char*) &tmpx, sizeof(Real));
				    ofs_mpi.write((char*) &tmpy, sizeof(Real));
				    ofs_mpi.write((char*) &tmpz, sizeof(Real));
				    ofs_mpi.write((char*) &tmpvol, sizeof(Real));
				    
				    for (int n = 0; n < comps.size(); ++n) {
				      ofs_mpi.write((char *) &(fab_array(x,y,z,comps[n])), sizeof(amrex::Real)); 
				    }
				    //					ofs_mpi.write(reinterpret_cast<const char*> (&fab_array(x,y,z,n), sizeof(amrex::Real));
				    //				  ofs_mpi<<fab_array(x,y,z,n) <<std::endl;
				}
		}
	      }
	    }
	    
	    //	    ofs<<fab<<std::endl;
	    //fab.writeOn(ofs_bin);
	    
        }
	
	//        if (lev>0)

    }

    }
    amrex::Finalize();
    return 0;
}

