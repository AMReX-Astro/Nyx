
#include "Nyx.H"
#include "Prob.H"

using namespace amrex;
using std::string;

void
Nyx::error_setup()
{
    std::string amr_prefix = "amr";
    ParmParse ppamr(amr_prefix);
    Vector<std::string> refinement_indicators;
    ppamr.queryarr("refinement_indicators",refinement_indicators,0,ppamr.countval("refinement_indicators"));

    if (refinement_indicators.size()==0)
        prob_errtags_default(errtags);

    else {
        for (int i=0; i<refinement_indicators.size(); ++i)
        {
        std::string ref_prefix = amr_prefix + "." + refinement_indicators[i];
        ParmParse ppr(ref_prefix);
        RealBox realbox;
        if (ppr.countval("in_box_lo")) {
          std::vector<Real> box_lo(AMREX_SPACEDIM), box_hi(AMREX_SPACEDIM);
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
          if(max_level >= 0)
            info.SetMaxLevel(max_level);
          else
            info.SetMaxLevel(0);
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
}

void
Nyx::manual_tags_placement (TagBoxArray&  /*tags*/,
                            const Vector<IntVect>& /*bf_lev*/)
{
}
