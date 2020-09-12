
#include "Nyx.H"
#include "Prob.H"
#include "Nyx_error_F.H"

using namespace amrex;
using std::string;

void
Nyx::error_setup()
{

#ifdef CXX_PROB
    std::string amr_prefix = "amr";
    ParmParse ppamr(amr_prefix);
    Vector<std::string> refinement_indicators;
    ppamr.queryarr("refinement_indicators",refinement_indicators,0,ppamr.countval("refinement_indicators"));
    if(refinement_indicators.size()==0)
        prob_errtags_default(errtags);
    else
    {
    for (int i=0; i<refinement_indicators.size(); ++i)
    {
    std::string ref_prefix = amr_prefix + "." + refinement_indicators[i];
	ParmParse ppr(ref_prefix);
	RealBox realbox;
	if (ppr.countval("in_box_lo")) {
	  std::vector<Real> box_lo(BL_SPACEDIM), box_hi(BL_SPACEDIM);
	  ppr.getarr("in_box_lo",box_lo,0,box_lo.size());
	  ppr.getarr("in_box_hi",box_hi,0,box_hi.size());
	  realbox = RealBox(&(box_lo[0]),&(box_hi[0]));
	}
	AMRErrorTagInfo info;
	if (realbox.ok()) {
	  info.SetRealBox(realbox);
	}
	if (ppr.countval("start_time") > 0) {
	  Real min_time; ppr.get("start_time",min_time);
	  info.SetMinTime(min_time);
	}
	if (ppr.countval("end_time") > 0) {
	  Real max_time; ppr.get("end_time",max_time);
	  info.SetMaxTime(max_time);
	}
	if (ppr.countval("max_level") > 0) {
	  int max_level; ppr.get("max_level",max_level);
	  info.SetMaxLevel(max_level);
	}
	if (ppr.countval("value_greater")) {
	  Real value; ppr.get("value_greater",value);
	  std::string field; ppr.get("field_name",field);
	  errtags.push_back(AMRErrorTag(value,AMRErrorTag::GREATER,field,info));
	}
	else if (ppr.countval("value_less")) {
	  Real value; ppr.get("value_less",value);
	  std::string field; ppr.get("field_name",field);
	  errtags.push_back(AMRErrorTag(value,AMRErrorTag::LESS,field,info));
	}
	else if (ppr.countval("vorticity_greater")) {
	  Real value; ppr.get("vorticity_greater",value);
	  const std::string field="mag_vort";
	  errtags.push_back(AMRErrorTag(value,AMRErrorTag::VORT,field,info));
	}
	else if (ppr.countval("adjacent_difference_greater")) {
	  Real value; ppr.get("adjacent_difference_greater",value);
	  std::string field; ppr.get("field_name",field);
	  errtags.push_back(AMRErrorTag(value,AMRErrorTag::GRAD,field,info));
	}
	else if (realbox.ok())
	  {
	    errtags.push_back(AMRErrorTag(info));
	  }
	else {
	  Abort(std::string("Unrecognized refinement indicator for " + refinement_indicators[i]).c_str());
	}
	}
	}

#else
    // The lines below define routines to be called to tag cells for error
    // estimation -- the arguments of each "add" call are:
    //   1. Name of variable (state variable or derived quantity) which will be
    //      passed into the Fortran subroutine.
    //   2. Number of ghost cells each array needs in each call to the Fortran
    //      subroutine.
    //   3. Type of Fortran subroutine -- this determines the argument list of
    //      the Fortran subroutine.  These types are pre-defined and are
    //      currently restricted to ErrorRec::Standard and ErrorRec::UseAverage.
    //   4. Name of Fortran subroutine.

    // The routine LAPLAC_ERROR uses the special evaluation of the second
    // derivative and can be called with any variable. Note that two ghost cells
    // are needed.

    //err_list.add("density", 2, ErrorRec::Standard,
    //             BL_FORT_PROC_CALL(TAG_LAPLAC_ERROR, tag_laplac_error));

    //err_list.add("pressure", 2, ErrorRec::Standard,
    //             BL_FORT_PROC_CALL(TAG_LAPLAC_ERROR, tag_laplac_error));

    //err_list.add("density", 1, ErrorRec::Standard,
    //             BL_FORT_PROC_CALL(TAG_DENERROR, tag_denerror));

    //err_list.add("Temp", 1, ErrorRec::Standard,
    //             BL_FORT_PROC_CALL(TAG_TEMPERROR, tag_temperror));

    //err_list.add("pressure", 1, ErrorRec::Standard,
    //             BL_FORT_PROC_CALL(TAG_PRESSERROR, tag_presserror));

    //err_list.add("x_velocity", 1, ErrorRec::Standard,
    //             BL_FORT_PROC_CALL(TAG_VELERROR,tag_velerror));

    //err_list.add("y_velocity", 1, ErrorRec::Standard,
    //             BL_FORT_PROC_CALL(TAG_VELERROR, tag_velerror));

    //err_list.add("z_velocity", 1, ErrorRec::Standard,
    //             BL_FORT_PROC_CALL(TAG_VELERROR, tag_velerror));

    //err_list.add("entropy", 1, ErrorRec::Standard,
    //             BL_FORT_PROC_CALL(TAG_ENTERROR, tag_enterror));

    err_list.add("total_density",1,ErrorRec::UseAverage,
                 BL_FORT_PROC_CALL(TAG_OVERDENSITY, tag_overdensity));

#endif
}

void
Nyx::manual_tags_placement (TagBoxArray&    tags,
                            const Vector<IntVect>& bf_lev)
{
}

