/*
 *  bundles.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 9/6/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#include <list>
#include <map>
#include <numeric>
#include <boost/math/distributions/binomial.hpp>

#include "common.h"
#include "bundles.h"
#include "scaffolds.h"

bool rc_stage = false;

using namespace std;
using boost::math::binomial;

//struct ScaffoldSorter
//{
//	ScaffoldSorter(RefSequenceTable& _rt) : rt(_rt) {} 
//	bool operator()(shared_ptr<Scaffold const> lhs, shared_ptr<Scaffold const> rhs)
//	{
//        assert (lhs);
//        assert (rhs);
//		const char* lhs_name = rt.get_name(lhs->ref_id());
//		const char* rhs_name = rt.get_name(rhs->ref_id());
//		int c = strcmp(lhs_name, rhs_name);
//		if (c != 0)
//		{
//			return c < 0;
//		}
//		if (lhs->left() != rhs->left())
//		{
//			return lhs->left() < rhs->left();
//		}
//        return false;
//	}
//	
//	RefSequenceTable& rt;
//};

struct ScaffoldSorter
{
	ScaffoldSorter(RefSequenceTable& _rt) : rt(_rt) {} 
	bool operator()(shared_ptr<Scaffold const> lhs, shared_ptr<Scaffold const> rhs)
	{
        //assert (lhs);
        //assert (rhs);
        if (!lhs || !rhs)
            return false;
		int lhs_order = rt.observation_order(lhs->ref_id());
        assert (lhs_order != -1);
		int rhs_order = rt.observation_order(rhs->ref_id());
        assert (rhs_order != -1);

		if (lhs_order != rhs_order)
		{
			return lhs_order < rhs_order;
		}
		if (lhs->left() != rhs->left())
		{
			return lhs->left() < rhs->left();
		}

		//make the sort stable by checking other stuff (needed for rc_data raw counts t_id)
		if (lhs->right() != rhs->right())
		{
			return (lhs->right() < rhs->right());
		}
		if (lhs->strand() != lhs->strand())
		{
			return (lhs->strand() < rhs->strand());
		}
		if (lhs->is_ref() && rhs->is_ref() && !lhs->annotated_trans_id().empty())
		{
			return (lhs->annotated_trans_id() < rhs->annotated_trans_id());
		}

        return false;
	}
	
	RefSequenceTable& rt;
};


int rc_cov_inc(int i) {
  return ++i;
}

void rc_update_tdata(HitBundle& bundle, const Scaffold& scaff,
	                          double cov, double fpkm) {
  if (!bundle.rcdata()) return;
  RC_BundleData& rc = *bundle.rcdata();
  if (rc.exons.size()==0) return;
  RC_ScaffData q(&scaff);
  set<RC_ScaffData>::iterator tdata = rc.tdata.find(q);
  if (tdata==rc.tdata.end()) {
	fprintf(stderr, "Error: cannot locate bundle ref. transcript %s (%d-%d)!\n",
		q.t_name.c_str(), q.l, q.r);
	return;
  }
  (*tdata).cov=cov;
  (*tdata).fpkm=fpkm;
}

FILE* rc_fwopen(const char* fname) {
 if (strcmp(fname,"-")==0) return stdout;
 string fpath(output_dir);
 fpath += "/";
 fpath += fname;
 fpath += ".ctab";
 FILE* fh=fopen(fpath.c_str(), "w");
 if (fh==NULL) {
   fprintf(stderr, "Error: cannot open file %s\n",
					fpath.c_str());
   exit(1);
   }
 return fh;
}


void rc_write_f2t(FILE* fh, map<uint, set<uint> >& f2t) {
  for (map<uint, set<uint> >::iterator m=f2t.begin(); m!=f2t.end(); ++m) {
    uint f_id=(*m).first;
    set<uint>& tset = (*m).second;
    for (set<uint>::iterator it=tset.begin();it!=tset.end();++it) {
 	 uint t_id = *it;
 	 fprintf(fh, "%u\t%u\n", f_id, t_id);
    }
  }
  fflush(fh);
}

void rc_write_fc(FILE* fh, const char* ref_name, RC_BundleData& rc, set<RC_Feature>& feats) {
  if (&rc.exons == &feats) {
	//writing counts for all exons in bundle
	for (set<RC_Feature>::iterator f=feats.begin(); f!=feats.end(); ++f) {
      const RC_Feature& exon = *f;
      //assert( exon.l >= rc.lmin );
      int L=exon.l-rc.lmin;
      int xlen=exon.r-exon.l;
      if (exon.l < rc.lmin) {
    	  //shouldn't be here
    	  if (exon.r<rc.lmin) continue;
    	  xlen-=(rc.lmin-exon.l);
    	  L=0;
      }
      if (rc.rmax<exon.r) {
    	  if (exon.l>rc.rmax) continue; //should never happen
    	  xlen-=(exon.r-rc.rmax);
      }
      int R=L+xlen;
      vector<int>::iterator xcov_begin;
      vector<int>::iterator xcov_end;
      vector<float>::iterator xmcov_begin;
      vector<float>::iterator xmcov_end;
      if (exon.strand=='+' || exon.strand=='.') {
         xcov_begin  = rc.f_cov.begin()+L;
         xcov_end = rc.f_cov.begin()+R;
         xmcov_begin = rc.f_mcov.begin()+L;
         xmcov_end = rc.f_mcov.begin()+R;
      } else {
        xcov_begin  = rc.r_cov.begin()+L;
        xcov_end = rc.r_cov.begin()+R;
        xmcov_begin = rc.r_mcov.begin()+L;
        xmcov_end = rc.r_mcov.begin()+R;
      }

      double avg = (double)accumulate(xcov_begin, xcov_end, 0) / xlen;
      vector<double> diff(xlen);
      transform(xcov_begin, xcov_end, diff.begin(),
                     bind2nd( minus<double>(), avg));
      double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
      double stdev = sqrt(sq_sum / xlen);

      double mavg = (double)accumulate(xmcov_begin, xmcov_end, 0) / xlen;
      vector<double> mdiff(xlen);
      transform(xmcov_begin, xmcov_end, mdiff.begin(),
                     bind2nd( minus<double>(), mavg));
      sq_sum = inner_product(mdiff.begin(), mdiff.end(), mdiff.begin(), 0.0);
      double mstdev = sqrt(sq_sum / xlen);
	  fprintf(fh,"%u\t%s\t%c\t%d\t%d\t%d\t%d\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\n",
		  exon.id, ref_name, exon.strand, exon.l+1, exon.r, exon.rcount,
		  exon.ucount, exon.mrcount, avg, stdev, mavg, mstdev);
	}
  } else { //introns
	for (set<RC_Feature>::iterator f=feats.begin(); f!=feats.end(); ++f) {
	  fprintf(fh,"%u\t%s\t%c\t%d\t%d\t%d\t%d\t%.2f\n",(*f).id, ref_name,
		  (*f).strand, (*f).l+1, (*f).r, (*f).rcount, (*f).ucount, (*f).mrcount);
	}
  }

  fflush(fh);
}

void rc_write_counts(const char* refname, HitBundle& bundle) {
 if (!bundle.rcdata()) return;
 RC_BundleData& rc = *bundle.rcdata();
 if (rc.exons.size()==0) return;
 //File: t_data.ctab
 //t_id tname chr strand start end num_exons gene_id gene_name cufflinks_cov cufflinks_fpkm
 for ( set<RC_ScaffData>::iterator si = rc.tdata.begin(); si!=rc.tdata.end();++si) {
  const RC_ScaffData& sd=*si;
  fprintf(rc.ftdata, "%u\t%s\t%c\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%lg\t%lg\n",
	  sd.t_id, refname, sd.strand, sd.l+1, sd.r, sd.t_name.c_str(),
	  sd.num_exons, sd.eff_len, sd.scaff->annotated_gene_id().c_str(),
	  sd.scaff->annotated_gene_name().c_str(), sd.cov, sd.fpkm);
 }//for each transcript
 fflush(rc.ftdata);

 //File: e_data.ctab
 //e_id chr gstart gend rcount ucount mrcount
 rc_write_fc(rc.fedata, refname, rc, rc.exons);
 //File: i_data.ctab
 //i_id chr gstart gend rcount ucount mrcount
 rc_write_fc(rc.fidata, refname, rc, rc.introns);

 // feature-to-transcript link files
 rc_write_f2t(rc.fe2t,  rc.e2t);
 rc_write_f2t(rc.fi2t,  rc.i2t);
}

//FIXME: needs refactoring
void load_ref_rnas(FILE* ref_mRNA_file, 
				   RefSequenceTable& rt,
				   vector<shared_ptr<Scaffold> >& ref_mRNAs,
				   bool loadSeqs,
				   bool loadFPKM) 
{
	if (loadSeqs)
		ProgressBar p_bar("Loading reference annotation and sequence.",0);
    else
        ProgressBar p_bar("Loading reference annotation.",0);

	GList<GSeqData> ref_rnas;
	
	// If the RefSequenceTable already has entries, we will sort the GTF records
	// according to their observation order.  Otherwise, we will sort the 
	// RefSequenceTable's records lexicographically.
	bool reorder_GTF_recs_lexicographically = false;
	if (rt.size() == 0)
	{
	  reorder_GTF_recs_lexicographically = true;
	}

	if (ref_mRNA_file)
	{
		gtf_tracking_verbose=cuff_verbose;
		read_transcripts(ref_mRNA_file, ref_rnas, true);
	}
	
	int last_gseq_id = -1;
	GFaSeqGet* faseq = NULL;
	GFastaHandler gfasta(fasta_dir.c_str());
	// Geo groups them by chr.
	if (ref_rnas.Count()>0) //if any ref data was loaded
	{
		for (int j = 0; j < ref_rnas.Count(); ++j) 
		{    //ref data is grouped by genomic sequence
			//const char* name = ref_rnas[j]->gseq_name;
			
			int f = 0;
			int r = 0;
			int u  = 0;
			GffObj* rna_p;
			RefID ref_id = rt.get_id(ref_rnas[j]->gseq_name, NULL);
			int f_count = ref_rnas[j]->mrnas_f.Count();
			int r_count = ref_rnas[j]->mrnas_r.Count();
			int u_count = ref_rnas[j]->umrnas.Count();
			
			while(!(f==f_count && r==r_count && u==u_count))
			{	
				CuffStrand strand;
				
				if (f < f_count)
				{
					rna_p = ref_rnas[j]->mrnas_f[f++];
					strand = CUFF_FWD;
				}
				else if (r < r_count) 
				{
					rna_p = ref_rnas[j]->mrnas_r[r++];
					strand = CUFF_REV;
				}
				else 
				{
					rna_p = ref_rnas[j]->umrnas[u++];
					strand = CUFF_STRAND_UNKNOWN;
				}

				GffObj& rna = *rna_p;

				if (loadSeqs && rna.gseq_id != last_gseq_id) //next chromosome
				{
					delete faseq;
					faseq = NULL;
					last_gseq_id = rna.gseq_id;
					faseq = gfasta.fetch(last_gseq_id);
					if (faseq==NULL)
					{
						fprintf(stderr,"This contig will not be bias corrected.\n");
					}
				}

				vector<AugmentedCuffOp> ops;
				for (int e = 0; e < rna.exons.Count(); ++e)
				{
					GffExon& ex = *(rna.exons[e]);
					ops.push_back(AugmentedCuffOp(CUFF_MATCH, ex.start - 1, ex.end - ex.start + 1));
					
					if (e + 1 < rna.exons.Count())
					{
						GffExon& next_ex = *(rna.exons[e+1]);
						ops.push_back(AugmentedCuffOp(CUFF_INTRON, ex.end, next_ex.start - ex.end - 1));
					}
				}
				
				Scaffold ref_scaff(ref_id, strand, ops, true);
				
				char* rna_seq = 0;
				int seqlen=0;
				if (loadSeqs && faseq){ 
					rna_seq = rna.getSpliced(faseq, false, &seqlen);
				}

				if (rna.getID())
					ref_scaff.annotated_trans_id(rna.getID());
				
				
				if (rna.getGeneID())
					ref_scaff.annotated_gene_id(rna.getGeneID());
				
				if (rna.getGeneName())
					ref_scaff.annotated_gene_name(rna.getGeneName());
				
				
				char* nearest_ref_match = rna.getAttr("nearest_ref");
				char* class_code = rna.getAttr("class_code");
				
				if (nearest_ref_match && class_code)
				{
					ref_scaff.nearest_ref_id(nearest_ref_match);
					ref_scaff.nearest_ref_classcode(*class_code);
				}
				
				char* protein_id = rna.getAttr("p_id");
				if (protein_id)
					ref_scaff.annotated_protein_id(protein_id);
				
				
				char* tss_id = rna.getAttr("tss_id");
				if (tss_id)
					ref_scaff.annotated_tss_id(tss_id);
				
				
				if (loadFPKM)
				{
					const char* expr = rna.getAttr("FPKM");
					if (expr!=NULL) {
						if (expr[0]=='"') expr++;
						ref_scaff.fpkm(strtod(expr, NULL));
					}
				}
				
				if (loadSeqs)
                {
					string rs = (rna_seq) ? rna_seq:"";
					std::transform(rs.begin(), rs.end(), rs.begin(), (int (*)(int))std::toupper);
					ref_scaff.seq(rs);
					GFREE(rna_seq);
				}
                
				shared_ptr<Scaffold> scaff(new Scaffold());
                *scaff = ref_scaff;
                assert (scaff);
				ref_mRNAs.push_back(scaff); 
			}
		}
        
        BOOST_FOREACH (shared_ptr<Scaffold> s, ref_mRNAs)
        {
            assert (s);
        }
        
        if (reorder_GTF_recs_lexicographically)
        {
            rt.order_recs_lexicographically();
        }
        
		ScaffoldSorter sorter(rt);
		sort(ref_mRNAs.begin(), ref_mRNAs.end(), sorter);
      if (emit_raw_counts) {
        uint cur_tid=0;
        uint cur_exon_id=0;
        uint cur_intron_id=0;
        set<RC_ScaffSeg> exons;
        set<RC_ScaffSeg> introns;
        //assign unique transcript IDs based on the sorted order
        RefID last_refid=0;
        for (vector<shared_ptr<Scaffold> >::iterator sit=ref_mRNAs.begin(); sit!=ref_mRNAs.end(); ++sit) {
          (*sit)->rc_init(++cur_tid);
          if ((*sit)->ref_id()!=last_refid) {
             exons.clear();
             introns.clear();
             last_refid=(*sit)->ref_id();
             }
        //debugOnly:
        //printf("%d\t%d\t%d\t%s\n",(*sit)->rc_tid(), (*sit)->left(), (*sit)->right(), 
        //  (*sit)->annotated_trans_id().c_str());
        //--^

          (*sit)->rc_addFeatures(cur_exon_id, exons, cur_intron_id, introns);
        }
        //debugOnly:
        //exit(1);
      }
	}
	delete faseq;
}


int HitBundle::_next_id = 0;

bool HitBundle::add_hit(const MateHit& hit)
{
	if (_final)
    {
		return false;
    }
	
	// Update the bounds on the span
	if (hit.left() < _leftmost)
		_leftmost = hit.left();
	if (hit.right() > _rightmost)
		_rightmost = hit.right();
	
	
	_hits.push_back(hit);
	return true;
}

struct HitlessScaffold
{
	bool operator()(shared_ptr<Scaffold> x)
	{
		return x->mate_hits().empty();
	}
};

bool unmapped_hit(const MateHit& x)
{
	return !(x.is_mapped());
}

void HitBundle::rc_store_t(Scaffold& scaff) {
 //if (this->ref_scaffolds().size()==0) return;
  if (!rc_stage) return;
  if (rc_data==NULL) {
	rc_init();
  }
  rc_data->addScaffold(scaff);
 //check this read alignment against ref exons and introns
}

struct COvlSorter {
  bool operator() (pair<int, const RC_Feature*> i,
      pair<int, const RC_Feature*> j) {
    return (i.first>j.first); //sort in decreasing order of overlap length
  }
} OvlSorter;


void rc_updateExonCounts(const RC_Feature* exon, const ReadHit* bh) {
  exon->rcount++;
  exon->mrcount += (bh->num_hits() > 1) ? (1.0/bh->num_hits()) : 1;
  if (bh->num_hits()<=1)  exon->ucount++;
}

void HitBundle::rc_count_hit(const ReadHit* bh) {
 if (rc_data==NULL) return;
 if (rc_data->tdata.size()==0) return;

 //check this read alignment against ref exons and introns
 char strand=(bh->source_strand()==CUFF_STRAND_UNKNOWN || bh->source_strand()==CUFF_BOTH) ?
         '.' : ((bh->source_strand()==CUFF_FWD) ? '+' : '-');
 int gstart=bh->left(); //alignment start position on the genome
 if (rc_data->f_cov.size()==0 && gstart<rc_data->lmin) {
   fprintf(stderr, "Warning: adjusting lmin coverage bundle from %d to %d !\n", int(rc_data->lmin), (int)gstart);
   rc_data->lmin=gstart;
 }
 if (bh->right()<rc_data->lmin) {
	 return; //hit outside coverage area
 }
 int gpos=bh->left(); //current genomic position
 int rlen=0; //read length, obtained here from the cigar string
 int segstart=gstart;
 vector<RC_Seg> rsegs;
 vector<RC_Seg> rintrons;
 for (size_t i = 0; i < bh->cigar().size(); ++i)
  {
   const CigarOp& op = bh->cigar()[i];
   switch (op.opcode) {
	  case MATCH:
        rlen += op.length;
        rc_data->updateCov(strand, bh->num_hits(), gpos, op.length);
      case DEL:
        gpos += op.length;
        break;
      case SOFT_CLIP:
      case INS:
        rlen += op.length;
        break;
      case REF_SKIP:
        //intron starts here
        //last exon ends here
        rsegs.push_back(RC_Seg(segstart, gpos) );
        int istart=gpos;
        gpos += op.length;
        segstart = gpos;
        rintrons.push_back(RC_Seg(istart, segstart));
        break;
   }
  } //foreach cigar op
  rsegs.push_back(RC_Seg(segstart, gpos));
  //now check rexons and rintrons with findExons() and findIntron()
  for (size_t i=0;i<rintrons.size();++i) {
   RC_FeatIt ri=rc_data->findIntron(rintrons[i].l, rintrons[i].r, strand);
   if (ri!=rc_data->introns.end()) {
   (*ri).rcount++;
   (*ri).mrcount += (bh->num_hits() > 1) ? (1.0/bh->num_hits()) : 1;
   if (bh->num_hits()==1)  (*ri).ucount++;
   }
  } //for each intron

   for (size_t i=0;i<rsegs.size();++i) {
     RC_FeatPtrSet ovlex=rc_data->findExons(rsegs[i].l, rsegs[i].r, strand);
     if (ovlex.size()==0) continue;
     if (ovlex.size()>1) {
       vector< pair<int, const RC_Feature*> > xovl; //overlapped exons sorted by decreasing overlap length
       for (RC_FeatPtrSet::iterator ox=ovlex.begin();ox != ovlex.end(); ++ox) {
         int ovlen=(*ox)->ovlen(rsegs[i].l, rsegs[i].r);
         if (ovlen>=5)
           xovl.push_back(pair<int, const RC_Feature*>(ovlen, *ox));
       }
       if (xovl.size()>1) {
          sort(xovl.begin(), xovl.end(), OvlSorter); //larger overlaps first
          //update the counts only for ref exons with max overlap to this segment
          int max_ovl=xovl.begin()->first;
          for (vector<pair<int, const RC_Feature*> >::iterator xo=xovl.begin();xo!=xovl.end();++xo) {
        	if (max_ovl - xo->first > 5 ) break; //more than +5 bases coverage for the other exons
        	rc_updateExonCounts(xo->second, bh);
          }
       } else if (xovl.size() == 1) {
         rc_updateExonCounts(xovl.begin()->second, bh);
       }
     } else {
       // 1 exon overlap only
       int ovlen=(*ovlex.begin())->ovlen(rsegs[i].l, rsegs[i].r);
       if (ovlen>=5) rc_updateExonCounts(*ovlex.begin(), bh);
     }

   } //for each read "exon"

}


bool HitBundle::add_open_hit(shared_ptr<ReadGroupProperties const> rg_props,
                             const ReadHit* bh,
							 bool expand_by_partner)
{
    assert (bh != NULL);
    
	_leftmost = min(_leftmost, bh->left());
	_ref_id = bh->ref_id();
    
	if (bh->is_singleton() || no_read_pairs)
	{
		_rightmost = max(_rightmost, bh->right());
		MateHit m(rg_props, bh->ref_id(), bh, NULL);
        if (m.right() - m.left() > max_gene_length)
        {
            fprintf(stderr, "Warning: hit is longer than max_gene_length, skipping\n");
            return false;
        }
		add_hit(m);
	}
	else
	{
        if (abs(bh->right() - bh->partner_pos()+1) > max_gene_length)
        {
            fprintf(stderr, "Warning: hit is longer than max_gene_length, skipping\n");
            return false;
        }
		if (expand_by_partner)
			_rightmost = max(max(_rightmost, bh->right()), bh->partner_pos()+1);
		OpenMates::iterator mi = _open_mates.find(bh->left());
		
		// Does this read hit close an open mate?
		if (mi == _open_mates.end())
		{
			// No, so add it to the list of open mates, unless we would
			// already have seen it's partner
			if(bh->left() <= bh->partner_pos())
			{
				MateHit open_hit(rg_props,
                                 bh->ref_id(), 
                                 bh, 
                                 NULL);
				
				pair<OpenMates::iterator, bool> ret;
				ret = _open_mates.insert(make_pair(bh->partner_pos(), 
												  list<MateHit>()));
				
				ret.first->second.push_back(open_hit);
			}
			else
			{
                // This should never happen during hit_driven or ref_guided bundling, and in the case of
                // ref_driven, this read clearly shouldn't map to any of the transcripts anyways.
                // Adding this hit would cause problems with multi-reads that straddle boundaries after assembly.
				// add_hit(MateHit(rg_props,bh->ref_id(), bh, NULL));
                return false;
			}
		}
		else
		{
			
			bool found_partner = false;
			// Maybe, see if we can find an ID match in the list of
			// open mates expecting a partner at this position
			for (list<MateHit>::iterator pi = mi->second.begin();
				 pi != mi->second.end();
				 ++pi)
			{
				MateHit& pm = *pi;
				
				if (pm.insert_id() == bh->insert_id())
				{
					// Found a partner?
					
					Scaffold L(MateHit(rg_props, bh->ref_id(), pm.left_alignment(), NULL));
					Scaffold R(MateHit(rg_props, bh->ref_id(), bh, NULL));
					
					bool strand_agree = L.strand() == CUFF_STRAND_UNKNOWN ||
					R.strand() == CUFF_STRAND_UNKNOWN ||
					L.strand() == R.strand();
					
					//bool orientation_agree = pm.left_alignment()->antisense_align() != bh->antisense_align();
					
					if (strand_agree && 
                        (!Scaffold::overlap_in_genome(L, R, olap_radius) ||
                         Scaffold::compatible(L,R)))
					{					
						pm.right_alignment(bh);
						add_hit(pm);
						mi->second.erase(pi);
						if (mi->second.empty())
							_open_mates.erase(mi);
						
						found_partner = true;
						break;
					}
				}
			}
			
			if (!found_partner)
			{
				// If we got here, couldn't actually close any mates with
				// this read hit, so open a new one, unless we can never
				// close this one
				if(bh->left() <= bh->partner_pos())
				{
					MateHit open_hit(rg_props, bh->ref_id(), bh, NULL);
					
					pair<OpenMates::iterator, bool> ret;
					ret = _open_mates.insert(make_pair(bh->partner_pos(), 
													  list<MateHit>()));
					
					ret.first->second.push_back(open_hit);
				}
				else
				{
                    // This should never happen during hit_driven or ref_guided bundling, and in the case of
                    // ref_driven, this read clearly shouldn't map to any of the transcripts anyways.
                    // Adding this hit would cause problems with multi-reads that straddle boundaries after assembly.
					// add_hit(MateHit(rg_props, bh->ref_id(), bh, NULL));
                    return false;
				}
			}
		}
	}
    return true;
}

void HitBundle::collapse_hits()
{
	::collapse_hits(_hits, _non_redundant);
    //_non_redundant = _hits;
}

void HitBundle::finalize_open_mates()
{
    // We don't want to split reads accross boundaries since this would only occur
    // in ref_driven mode and the read shouldn't map to any of the references in this case.

    for(OpenMates::iterator itr = _open_mates.begin(); itr != _open_mates.end(); ++itr)
    {
        BOOST_FOREACH (MateHit& hit,  itr->second)
        {
            delete hit.left_alignment();
            delete hit.right_alignment();
        }
    }
    _open_mates.clear();
}

void HitBundle::remove_hitless_scaffolds()
{
	vector<shared_ptr<Scaffold> >::iterator new_end = remove_if(_ref_scaffs.begin(),
												   _ref_scaffs.end(),
												   HitlessScaffold());
	_ref_scaffs.erase(new_end, _ref_scaffs.end());	
}



void HitBundle::combine(const vector<HitBundle*>& in_bundles,
                        HitBundle& out_bundle)
{
    out_bundle._hits.clear();
    out_bundle._non_redundant.clear();
    out_bundle._ref_scaffs.clear();
    
    for (size_t i = 1; i < in_bundles.size(); ++i)
    {
        assert(in_bundles[i]->ref_id() == in_bundles[i-1]->ref_id());
    }
    
    // Merge  hits
    vector<size_t> indices(in_bundles.size(),0);
    while(true)
    {
        int next_bundle = -1;
        const MateHit* next_hit=NULL; 
        for(size_t i = 0; i < in_bundles.size(); ++i)
        {
            const vector<MateHit>& curr_hits = in_bundles[i]->hits();
            
            if (indices[i] == curr_hits.size())
                continue;
            
            const MateHit* curr_hit = &curr_hits[indices[i]];
            
            if (next_bundle == -1 || mate_hit_lt(*curr_hit, *next_hit))
            {
                next_bundle = i;
                next_hit = curr_hit;
            }
        }
        
        if(next_bundle==-1)
            break;
        
        out_bundle._hits.push_back(*next_hit);
        indices[next_bundle]++;
    }
    
    // Merge collapsed hits
    indices = vector<size_t>(in_bundles.size(), 0);
    while(true)
    {
        int next_bundle = -1;
        const MateHit* next_hit = NULL; 
        for(size_t i = 0; i < in_bundles.size(); ++i)
        {
            const vector<MateHit>& curr_non_redundant_hits = in_bundles[i]->non_redundant_hits();
            
            if (indices[i] == curr_non_redundant_hits.size())
                continue;
            
            const MateHit* curr_hit = &curr_non_redundant_hits[indices[i]];
            
            if (next_bundle == -1 || mate_hit_lt(*curr_hit, *next_hit))
            {
                next_bundle = i;
                next_hit = curr_hit;
            }
        }
        
        if(next_bundle==-1)
            break;
        
        out_bundle._non_redundant.push_back(*next_hit);
        indices[next_bundle]++;
    }
    
    for(size_t i = 0; i < in_bundles.size(); ++i)
    {
        for (size_t j = 0; j < in_bundles[i]->_ref_scaffs.size(); ++j)
        {
            in_bundles[i]->_ref_scaffs[j]->clear_hits();
        }
    }
    
    // Merge ref scaffolds
    indices = vector<size_t>(in_bundles.size(), 0);
    while(true)
    {
        int next_bundle = -1;
        shared_ptr<Scaffold> next_scaff; 
        for(size_t i = 0; i < in_bundles.size(); ++i)
        {
            const vector<shared_ptr<Scaffold> >& curr_scaffs = in_bundles[i]->_ref_scaffs;
            
            if (indices[i] == curr_scaffs.size())
                continue;
            
            shared_ptr<Scaffold> curr_scaff = curr_scaffs[indices[i]];
            
            if (next_bundle == -1 || scaff_lt_rt_oplt(*curr_scaff, *next_scaff))
            {
                next_bundle = i;
                next_scaff = curr_scaff;
            }
        }
        
        if(next_bundle==-1)
            break;
        
        if (out_bundle._ref_scaffs.size()==0 || out_bundle._ref_scaffs.back()->annotated_trans_id() != next_scaff->annotated_trans_id()) 
            out_bundle.add_ref_scaffold(next_scaff);
        indices[next_bundle]++;
    }
	
    out_bundle.finalize(true); // true means everything is already sorted, etc.
    out_bundle._num_replicates = (int)in_bundles.size();
}


void HitBundle::finalize(bool is_combined)
{
	_final = true;
    if (!is_combined)
	{
        // only perform read skipping on primary bundles 
        // (i.e. don't do it on bundles we're making by combining two or more other bundles)
        size_t num_skipped = _hits.size() * read_skip_fraction;
        if (num_skipped > 0 && num_skipped < _hits.size())
        {
            random_shuffle(_hits.begin(), _hits.end());
            for (int i = (int)_hits.size() - num_skipped; i >= 0 && i < (int)_hits.size(); ++i)
            {
                delete _hits[i].left_alignment();
                _hits[i].left_alignment(NULL);
                
                delete _hits[i].right_alignment();
                _hits[i].right_alignment(NULL);
            }
            _hits.resize(_hits.size() - num_skipped);
            is_combined = false;
        }
        else if (num_skipped >= _hits.size())
        {
            for (size_t i = 0; i < _hits.size(); ++i)
            {
                delete _hits[i].left_alignment();
                delete _hits[i].right_alignment();
            }
            _hits.clear();
        }

		sort(_hits.begin(), _hits.end(), mate_hit_lt);
        if (cond_prob_collapse)
        {
            collapse_hits();
        }
        else
        {
            BOOST_FOREACH (MateHit& hit, _hits)
            {
                hit.incr_collapse_mass(hit.internal_scale_mass());
            }
            _non_redundant = _hits;
            
        }
		sort(_ref_scaffs.begin(), _ref_scaffs.end(), scaff_lt_rt_oplt_sp);
		vector<shared_ptr<Scaffold> >::iterator new_end = unique(_ref_scaffs.begin(), 
												_ref_scaffs.end(),
												StructurallyEqualScaffolds());
		_ref_scaffs.erase(new_end, _ref_scaffs.end());
        vector<shared_ptr<Scaffold> >(_ref_scaffs).swap(_ref_scaffs);
	}
	
    for (size_t j = 0; j < _ref_scaffs.size(); ++j)
	{
		_ref_scaffs[j]->clear_hits();
	}
    
    _compatible_mass = 0.0;
    
	for (size_t i = 0; i < _hits.size(); ++i)
	{
		MateHit& hit = _hits[i];
		
		Scaffold hs(hit);
		
        if (i >= 1)
        {
            assert (hit.ref_id() == _hits[i-1].ref_id());
        }
		hit.is_mapped(false);
		for (size_t j = 0; j < _ref_scaffs.size(); ++j)
		{
			// add hit only adds if the hit is structurally compatible
			if (_ref_scaffs[j]->contains(hs))
			{
				bool added = _ref_scaffs[j]->add_hit(&hit);
                if (added)
                    hit.is_mapped(true);
			}
		}
        if (hit.is_mapped())
        {
            _compatible_mass += hit.internal_scale_mass();
        }
	}
    
}

void print_sort_error(const char* last_chr_name, 
                      int last_chr_pos, 
                      const char* bh_name, 
                      int bh_pos)
{
    fprintf(stderr, "\nError: this SAM file doesn't appear to be correctly sorted!\n");
    fprintf(stderr, "\tcurrent hit is at %s:%d, last one was at %s:%d\n", 
            bh_name,
            bh_pos,
            last_chr_name,
            last_chr_pos);
    fprintf(stderr, "Cufflinks requires that if your file has SQ records in\nthe SAM header that they appear in the same order as the chromosomes names \nin the alignments.\nIf there are no SQ records in the header, or if the header is missing,\nthe alignments must be sorted lexicographically by chromsome\nname and by position.\n \n");
}


double BundleFactory::next_valid_alignment(const ReadHit*& bh)
{
    const char* hit_buf;
	size_t hit_buf_size = 0;
    bh = NULL;
    
	// Keep track of mass of hits we skip
	double raw_mass = 0; 
	
    while (true)
    {
    
        if (!_hit_fac->next_record(hit_buf, hit_buf_size))
            break;
        
        ReadHit tmp;
        if (!_hit_fac->get_hit_from_buf(hit_buf, tmp, false))
            continue;
        
		if (tmp.ref_id() == 12638153115695167477)  // corresponds to SAM "*" under FNV hash. unaligned read record 
            continue;
        
		raw_mass += tmp.mass();
		
        if (_hit_fac->ref_table().get_name(tmp.ref_id())==NULL) // unaligned read record (!?)
            continue;
            
        if (spans_bad_intron(tmp))
            continue;
        
        int order = _hit_fac->ref_table().observation_order(tmp.ref_id());
        if (_prev_pos != 0)
        {
            int prev_order = _hit_fac->ref_table().observation_order(_prev_ref_id);
            
            if (prev_order > order || (prev_order == order && _prev_pos > tmp.left()))
            {
                const char* bh_chr_name = _hit_fac->ref_table().get_name(tmp.ref_id());
                const char* last_bh_chr_name = _hit_fac->ref_table().get_name(_prev_ref_id);
                                
                print_sort_error(last_bh_chr_name, 
                                 _prev_pos, 
                                 bh_chr_name, 
                                 tmp.left());
                exit(1);
            }
        }
        
        _prev_ref_id = tmp.ref_id();
        _prev_pos = tmp.left();
        
        bool hit_within_mask = false;
        
        // We want to skip stuff that overlaps masked GTF records, so 
        // sync up the masking chromosome
        if (!mask_gtf_recs.empty() && 
            next_mask_scaff != mask_gtf_recs.end() &&
            (*next_mask_scaff)->ref_id() != tmp.ref_id())
        {
            bool found_scaff = false;
            vector<shared_ptr<Scaffold> >::iterator curr_mask_scaff = mask_gtf_recs.begin();
            for (size_t i = 0; i < _mask_scaff_offsets.size(); ++i)
            {
                if (_mask_scaff_offsets[i].first == tmp.ref_id())
                {
                    curr_mask_scaff = _mask_scaff_offsets[i].second;
                    found_scaff = true;
                    break;
                }
            }
            
            next_mask_scaff = curr_mask_scaff;
        }
        
        //check that we aren't sitting in the middle of a masked scaffold
        while (next_mask_scaff != mask_gtf_recs.end() && 
               (*next_mask_scaff)->ref_id() == tmp.ref_id() &&
               (*next_mask_scaff)->right() <= tmp.left())
        {
            if ((*next_mask_scaff)->left() >= tmp.left())
            {
                break;
            }
            
            next_mask_scaff++;
        }
        
        if (next_mask_scaff != mask_gtf_recs.end() &&
            (*next_mask_scaff)->ref_id() == tmp.ref_id() &&
            (*next_mask_scaff)->left() <= tmp.left() &&
            (*next_mask_scaff)->right() >= tmp.right())
        {
            hit_within_mask = true;
        }
        
        if (hit_within_mask)
            continue;
        
        // if the user's asked for read trimming, do it here.
        if (trim_read_length > 0)
        {
            tmp.trim(trim_read_length);
        }
        
        bh = new ReadHit(tmp);
        
        break;
    }
    
    return raw_mass;
}

double BundleFactory::rewind_hit(const ReadHit* rh)
{
	double mass = rh->mass();
	delete rh;
	_hit_fac->undo_hit();
	return mass;
}

bool BundleFactory::next_bundle_hit_driven(HitBundle& bundle)
{
	const ReadHit* bh = NULL;
    
    bool skip_read = false;
    
	while(bh == NULL)
	{
		if (!_hit_fac->records_remain())
		{
			return false;
		}
        
        // If we are randomly throwing out reads, check to see
        // whether this one should be kept.
        if (bundle.hits().size() >= max_frags_per_bundle)
        {
            skip_read = true;
            next_valid_alignment(bh);
        }
        else
        {
            double raw_mass = next_valid_alignment(bh);
            if (bh && bh->num_hits() > max_frag_multihits)
            {
                skip_read = true;
            }
            else
            {
                bundle.add_raw_mass(raw_mass);
            }
        }
	}
	
	if ((skip_read || !bundle.add_open_hit(read_group_properties(), bh)) && bh != NULL)
    {
        delete bh;
        bh = NULL;
    }
	_expand_by_hits(bundle);

    assert(bundle.left() != -1);    
	bundle.finalize_open_mates();
	bundle.finalize();
    assert(bundle.right() != -1);
    
    return true;
}

bool BundleFactory::next_bundle_ref_driven(HitBundle& bundle)
{
	if (next_ref_scaff == ref_mRNAs.end())
	{
		const ReadHit* bh = NULL;
		while(_hit_fac->records_remain())
		{
            double raw_mass = next_valid_alignment(bh);
            if (bundle.hits().size() < max_frags_per_bundle)
            {
                if (bh && bh->num_hits() > max_frag_multihits)
                {
                    
                }
                else
                {
                    bundle.add_raw_mass(raw_mass);
                }
                if (bh) { delete bh; }
            }
            else
            {
                delete bh;
                bh = NULL;
            }
		    
		}
		bundle.finalize();
		return false;
	}
	
	bundle.add_ref_scaffold(*next_ref_scaff);
	++next_ref_scaff;
    
	_expand_by_refs(bundle);
	
	// The most recent RefID and position we've seen in the hit stream
	RefID last_hit_ref_id_seen = 0;
	int last_hit_pos_seen = 0;

	if (rc_stage) {
	  bundle.rc_finalize_refs();
	}

	// include hits that lay within the bundle interval
	while(true)
	{		
		const ReadHit* bh = NULL;
        
        bool skip_read = false;
		// If we are randomly throwing out reads, check to see
        // whether this one should be kept.
        if (bundle.hits().size() >= max_frags_per_bundle)
        {
            next_valid_alignment(bh);
            skip_read = true;
        }
        else
        {
            double raw_mass = next_valid_alignment(bh);
            if (bh && bh->num_hits() > max_frag_multihits)
            {
                skip_read = true;
            }
            else
            {
                bundle.add_raw_mass(raw_mass);
            }
        }
        
        if (bh == NULL)
        {
			if (_hit_fac->records_remain())
				continue;
			else
				break;
        }
        
		last_hit_ref_id_seen = bh->ref_id();
		last_hit_pos_seen = bh->left();
		
		// test if the hit stream needs to catch up or has gone too far based on ref_id
		if (bh->ref_id() != bundle.ref_id())
		{
			int bh_chr_order = _hit_fac->ref_table().observation_order(bh->ref_id());
			int bundle_chr_order = _hit_fac->ref_table().observation_order(bundle.ref_id());
			
			if (bh_chr_order < bundle_chr_order) // the hit stream has not caught up, skip
			{
				delete bh;
                bh = NULL;
				continue; 
			}
			else // the hit stream has gone too far, rewind and break
			{
                rewind_hit(bh);
                bh = NULL;
                break;
			}
		}
        
        if (bh == NULL) // the hit stream has gone too far, break
            break;
		
        if (bh->left() >= bundle.left() && bh->right() <= bundle.right())
		{
            if (emit_raw_counts && rc_stage) {
              bundle.rc_count_hit(bh);
            }
            if (skip_read)
            {
                delete bh;
                bh = NULL;
            }
			else
            {
                if (!bundle.add_open_hit(read_group_properties(), bh, false))
                {
                    delete bh;
                    bh = NULL;
                }
            }
		}
		else if (bh->left() >= bundle.right())
		{
            if (skip_read == false)
            {
                bundle.rem_raw_mass(rewind_hit(bh));
                bh = NULL;
            }
            else
            {
                delete bh;
                bh = NULL;
            }
			break;
		}
	    else
        {
            // It's not within the bundle bounds, but it's also not past the 
            // right end, so skip it.
            delete bh;
            bh = NULL;
        }
        
        if (skip_read == true && bh != NULL)
        {
            delete bh;
            bh = NULL;
        }
	}
	
    assert(bundle.left() != -1);
    bundle.finalize_open_mates();
	bundle.finalize();
    assert(bundle.right() != -1);
    
    return true;
}

// NOTE: does not support read skipping yet or max hits per bundle yet.
bool BundleFactory::next_bundle_ref_guided(HitBundle& bundle)
{
	
	if (next_ref_scaff == ref_mRNAs.end())
	{
		return next_bundle_hit_driven(bundle);
	}
	
	const ReadHit* bh = NULL;
	while(bh == NULL)
	{
		if (!_hit_fac->records_remain())
		{
			return next_bundle_ref_driven(bundle);
		}
		bundle.add_raw_mass(next_valid_alignment(bh));
	}
	
	if (bh->ref_id() != (*next_ref_scaff)->ref_id())
	{
		int bh_chr_order = _hit_fac->ref_table().observation_order(bh->ref_id());
		int scaff_chr_order = _hit_fac->ref_table().observation_order((*next_ref_scaff)->ref_id());
		
		bundle.rem_raw_mass(rewind_hit(bh));
		bh = NULL;
        
		if (bh_chr_order < scaff_chr_order)
		{
			return next_bundle_hit_driven(bundle);
		}
		else
		{
			return next_bundle_ref_driven(bundle);
		}
	}
		
	if (bh->left() < (*next_ref_scaff)->left())
	{
		if (!bundle.add_open_hit(read_group_properties(), bh))
        {
            delete bh;
            bh = NULL;
        }
	}
	else 
	{
		bundle.rem_raw_mass(rewind_hit(bh));
        bh = NULL;
        
		bundle.add_ref_scaffold(*next_ref_scaff);
		next_ref_scaff++;
		_expand_by_refs(bundle);
	}
	
	while(_expand_by_hits(bundle) || 
		  _expand_by_refs(bundle)) {}
	
	assert(bundle.left() != -1);    
	bundle.finalize_open_mates();
	bundle.finalize();
	assert(bundle.right() != -1);
	
	return true; 
}

// expand the bundle interval as far as needed to include the overlapping
// chain of reference transcripts that also overlap the initial bundle
// interval
bool BundleFactory::_expand_by_refs(HitBundle& bundle)
{
	int initial_right = bundle.right();
	while(next_ref_scaff < ref_mRNAs.end())
	{		
		assert(bundle.ref_id() != (*next_ref_scaff)->ref_id() || (*next_ref_scaff)->left() >= bundle.left());
		if (bundle.ref_id() == (*next_ref_scaff)->ref_id()
			&& overlap_in_genome((*next_ref_scaff)->left(),(*next_ref_scaff)->right(),bundle.left(), bundle.right()))
		{
			bundle.add_ref_scaffold(*next_ref_scaff++);
		}
		else 
		{
			break;
		}		
	}

	
	return (bundle.right() > initial_right);
}

// expand bundle by chaining overlapping hits
bool BundleFactory::_expand_by_hits(HitBundle& bundle)
{
	int initial_right = bundle.right();
	while(true)
	{
        bool skip_read = false;
        const ReadHit* bh = NULL;
        
        double raw_mass = next_valid_alignment(bh);
        if (bh && bh->num_hits() > max_frag_multihits)
        {
            skip_read = true;
        }
        else
        {
            bundle.add_raw_mass(raw_mass);
        }

		if (bh == NULL)
		{
			if (_hit_fac->records_remain())
			{
				continue;
			}
			else
			{
				break;
			}	
		}
		
		if (bh->ref_id() == bundle.ref_id() && bh->left() < bundle.right() + olap_radius)
		{			
			if (skip_read || !bundle.add_open_hit(read_group_properties(), bh))
            {
                delete bh;
                bh = NULL;
            }
		}
		else
		{
			bundle.rem_raw_mass(rewind_hit(bh));

			break;
		}
	}
	
	return (bundle.right() > initial_right);
}

bool BundleFactory::next_bundle(HitBundle& bundle)
{    
#if ENABLE_THREADS
    boost::mutex::scoped_lock lock(_factory_lock);
#endif
	switch(_bundle_mode)
	{
		case HIT_DRIVEN:
            _curr_bundle++;
			return next_bundle_hit_driven(bundle);
			break;
		case REF_DRIVEN:
            _curr_bundle++;
			return next_bundle_ref_driven(bundle);
			break;
		case REF_GUIDED:
            _curr_bundle++;
			return next_bundle_ref_guided(bundle);
			break;
	}
	return false;
}


struct IntronSpanCounter
{
	IntronSpanCounter() : left_reads(0), little_reads(0), total_reads(0), multimap_reads(0), fwd_strand_frags(0) {}
	size_t left_reads;
	size_t little_reads; // small span overhang
	size_t total_reads;
	size_t multimap_reads;
    size_t fwd_strand_frags;
	vector<size_t> hist;
};

typedef map<AugmentedCuffOp, IntronSpanCounter> IntronCountTable;

void count_introns_in_read(const ReadHit& read,
						   IntronCountTable& intron_counts)
{
	const vector<CigarOp>& cig = read.cigar();
	
	int read_len = read.read_len();
	int small_anchor = (int)floor(read_len * small_anchor_fraction);
	
	int r_left = 0;
	int g_left = read.left();
	
	for (size_t i = 0; i < cig.size(); ++i)
	{
		assert(cig[i].length >= 0);
		switch(cig[i].opcode)
		{
			case MATCH:
				//ops.push_back(AugmentedCuffOp(CUFF_MATCH, g_left, cig[i].length));
				g_left += cig[i].length;
				r_left += cig[i].length;
				break;
				
			case REF_SKIP:
			{	
				AugmentedCuffOp intron(CUFF_INTRON, g_left, cig[i].length);
				pair<IntronCountTable::iterator, bool> ins_itr;
				ins_itr = intron_counts.insert(make_pair(intron, IntronSpanCounter()));
				IntronCountTable::iterator itr = ins_itr.first;
				itr->second.total_reads++;
				
				if (read.num_hits() > 10)
				{
					itr->second.multimap_reads++;
				}
				
				if ( r_left <= small_anchor || (read_len - r_left) < small_anchor)
				{
					itr->second.little_reads++;
				}
                
                if (read.source_strand() == CUFF_FWD)
                {
                    //itr->second.fwd_strand_frags;
                }
                else 
                {
                    assert(read.source_strand() == CUFF_REV);
                }

				
				vector<size_t>& hist = itr->second.hist;
				if (hist.size() < (size_t)read_len)
				{
					size_t num_new_bins = read_len - hist.size();
					size_t new_left_bins = (size_t)floor(num_new_bins / 2.0);
					size_t new_right_bins = (size_t)ceil(num_new_bins / 2.0);
					hist.insert(hist.begin(), new_left_bins, 0);
					hist.insert(hist.end(), new_right_bins, 0);
				}
				
				assert (r_left < hist.size());
				hist[r_left]++;
				//ops.push_back(AugmentedCuffOp(CUFF_INTRON, g_left, cig[i].length));
				g_left += cig[i].length;
				break;
			}
				
			case SOFT_CLIP:
				g_left += cig[i].length;
				break;
            case HARD_CLIP:
				break;
            case INS:
                g_left -= cig[i].length;
                break;
            case DEL:
                g_left += cig[i].length;
                break;
			default:
				assert(false);
				break;
		}
	}
}

void minor_introns(int bundle_length,
				   int bundle_left,
				   const IntronCountTable& intron_counts,
				   vector<AugmentedCuffOp>& bad_introns,
				   double fraction)

{
	for(IntronCountTable::const_iterator itr = intron_counts.begin();
		itr != intron_counts.end(); 
		++itr)
	{
		pair<AugmentedCuffOp, IntronSpanCounter> itr_cnt_pair = *itr;
		const IntronSpanCounter itr_spans = itr_cnt_pair.second;
		
		double doc = itr_spans.total_reads;
		
		for (IntronCountTable::const_iterator itr2 = intron_counts.begin();
			 itr2 != intron_counts.end(); 
			 ++itr2)
		{	
			if (itr == itr2 ||
				!AugmentedCuffOp::overlap_in_genome(itr->first, itr2->first))
			{
				continue;
			}
			
			pair<AugmentedCuffOp, IntronSpanCounter> itr2_cnt_pair = *itr2;
			const IntronSpanCounter itr2_spans = itr2_cnt_pair.second;
			
			double thresh = itr2_spans.total_reads * fraction;
			if (doc < thresh)
			{
				//#if verbose_msg
				//							fprintf(stderr, "\t Filtering intron (due to overlap) %d - %d: %f thresh %f\n", itr->first.first, itr->first.second, doc, bundle_avg_thresh);
				//#endif	
				bool exists = binary_search(bad_introns.begin(), 
											bad_introns.end(), 
											itr->first);
				if (!exists)
				{
					verbose_msg("Filtering intron %d-%d spanned by %lu reads based on overlap with much more abundant intron: %d-%d spanned by %lu reads\n", 
							itr->first.g_left(), 
							itr->first.g_right(), 
							itr->second.total_reads,
							itr2->first.g_left(), 
							itr2->first.g_right(), 
							itr2->second.total_reads);
					
					bad_introns.push_back(itr->first);
					sort(bad_introns.begin(), bad_introns.end());
				}
			}
            
//            if ((itr->second.fwd_strand_frags == 0 &&
//                 itr2->second.fwd_strand_frags != 0) ||
//                (itr2->second.fwd_strand_frags == 0 &&
//                 itr->second.fwd_strand_frags != 0))
//            {
//                int itr1_L = itr->first.g_left();
//                int itr1_R = itr->first.g_right();
//                int itr2_L = itr2->first.g_left();
//                int itr2_R = itr2->first.g_right();
//                
//                if (abs(itr1_L - itr2_L) < 25 && abs(itr1_R - itr2_R) < 25)
//                {
//                    int a = 3;
//                }
//            }
		}
	}
}

void multimapping_introns(int bundle_length,
						  int bundle_left,
						  const IntronCountTable& intron_counts,
						  vector<AugmentedCuffOp>& bad_introns,
						  double fraction)

{
	for(IntronCountTable::const_iterator itr = intron_counts.begin();
		itr != intron_counts.end(); 
		++itr)
	{
		pair<AugmentedCuffOp, IntronSpanCounter> itr_cnt_pair = *itr;
		const IntronSpanCounter itr_spans = itr_cnt_pair.second;
		
		double doc = itr_spans.total_reads;
		double multi = itr_spans.multimap_reads;
		
		double multi_fraction = multi / doc;
		
		if (multi_fraction > fraction)
		{
			bool exists = binary_search(bad_introns.begin(), 
										bad_introns.end(), 
										itr->first);
			if (!exists)
			{
				verbose_msg("Filtering intron %d-%d spanned by %lu reads because %lg percent are multireads.\n", 
						itr->first.g_left(), 
						itr->first.g_right(), 
						itr->second.total_reads,
						multi_fraction * 100);
				
				bad_introns.push_back(itr->first);
				sort(bad_introns.begin(), bad_introns.end());
			}
		}
	}
}


void identify_bad_splices(const HitBundle& bundle, 
						  BadIntronTable& bad_splice_ops)
{
	// Tracks, for each intron, how many reads
	IntronCountTable intron_counts;
	
	RefID ref_id = bundle.ref_id();
	
	pair<BadIntronTable::iterator, bool> ins_itr;
	ins_itr = bad_splice_ops.insert(make_pair(ref_id, vector<AugmentedCuffOp>()));
	vector<AugmentedCuffOp>& bad_introns = ins_itr.first->second;
	
	BOOST_FOREACH (const MateHit& hit, bundle.hits())
	{
		if (hit.left_alignment())
		{
			count_introns_in_read(*hit.left_alignment(), intron_counts);
		}
		if (hit.right_alignment())
		{
			count_introns_in_read(*hit.right_alignment(), intron_counts);
		}
	}
	
	minor_introns(bundle.length(), bundle.left(), intron_counts, bad_introns, min_isoform_fraction);
	// [Geo] disable filtering of multi-mapped introns:
  // multimapping_introns(bundle.length(), bundle.left(), intron_counts, bad_introns, 0.5);
	for (IntronCountTable::iterator itr = intron_counts.begin();
		 itr != intron_counts.end();
		 ++itr)
	{
		if (binary_search(bad_introns.begin(), 
						  bad_introns.end(), 
						  itr->first))
		{
			continue;
		}
		pair<AugmentedCuffOp, IntronSpanCounter> cnt_pair = *itr;
		try
		{
			const IntronSpanCounter spans = cnt_pair.second;
			
			//			binomial read_half_dist(spans.total_reads, success_fraction);
			//			double left_side_p = cdf(read_half_dist, spans.total_reads - spans.left_reads);
			//			double right_side_p = cdf(complement(read_half_dist, spans.left_reads));
			
			
			double success = 2 * small_anchor_fraction;
			
			binomial read_half_dist(spans.total_reads, success);
			double right_side_p;
			
			// right_side_p describes the chance that we'd observe at least 
			// this many small overhang reads by chance with an unbiased 
			// distribution over a normal (e.g. non-artifact) junction
			if (spans.little_reads > 0)
			{
				right_side_p = 1.0 - cdf(read_half_dist, spans.little_reads - 1);
			}
			else 
			{
				right_side_p = 1.0;
			}
			
			double left_side_p = 0;
			double expected = success * spans.total_reads;

            //double excess = spans.little_reads - expected;
			
			// left_side_p describes the chance that we'd observe this few or
			// fewer small overhang reads by chance with an unbiased 
			// distribution over a normal (e.g. non-artifact) junction
			if (spans.little_reads > 0)
			{
				left_side_p = cdf(read_half_dist, spans.little_reads);
			}
			else 
			{
				left_side_p = cdf(read_half_dist, 0);
			}
			
			//double alpha = 0.05;
			//double right_side_p = 0;
			
			// Two-tailed binomial test:
//			if (left_side_p < (binomial_junc_filter_alpha / 2.0) || 
//				right_side_p < (binomial_junc_filter_alpha / 2.0))
			// One-tailed binomial test
			
			bool filtered = false;
			
			const IntronSpanCounter& counter = itr->second;
			
			if (right_side_p < (binomial_junc_filter_alpha))
			{
				double overhang_ratio = counter.little_reads / (double) counter.total_reads;
				if (counter.total_reads < 100 || overhang_ratio >= 0.50)
				{
					verbose_msg("Filtering intron %d-%d spanned by %lu reads (%lu low overhang, %lg expected) left P = %lg, right P = %lg\n", 
							itr->first.g_left(), 
							itr->first.g_right(), 
							itr->second.total_reads, 
							itr->second.little_reads, 
							expected,
							left_side_p,
							right_side_p);
					filtered = true;
					
					bool exists = binary_search(bad_introns.begin(), 
												bad_introns.end(), 
												itr->first);
					if (!exists)
					{
						bad_introns.push_back(itr->first);
						sort(bad_introns.begin(), bad_introns.end());
					}
				}
			}
			
			vector<size_t> hist = itr->second.hist;
			if (itr->second.total_reads > 1000)
			{
				sort(hist.begin(), hist.end());
				size_t median = (size_t)floor(hist.size() / 2);
				if (median <= hist.size() && hist[median] == 0)
				{
					verbose_msg("Filtering intron %d-%d spanned by %lu reads (%lu low overhang, %lg expected) left P = %lg, right P = %lg\n", 
							itr->first.g_left(), 
							itr->first.g_right(), 
							itr->second.total_reads, 
							itr->second.little_reads, 
							expected,
							left_side_p,
							right_side_p);
					
					filtered = true;
					
					bool exists = binary_search(bad_introns.begin(), 
												bad_introns.end(), 
												itr->first);
					if (!exists)
					{
						bad_introns.push_back(itr->first);
						sort(bad_introns.begin(), bad_introns.end());
					}
				}
			}
			
			if (!filtered)
			{
				verbose_msg("Accepting intron %d-%d spanned by %lu reads (%lu low overhang, %lg expected) left P = %lg, right P = %lg\n", 
						itr->first.g_left(), 
						itr->first.g_right(), 
						itr->second.total_reads, 
						itr->second.little_reads, 
						expected,
						left_side_p,
						right_side_p);
				
			}
		}
		
		
		catch(const std::exception& e)
		{
			//
			/*`
			 [#coinflip_eg_catch]
			 It is always essential to include try & catch blocks because
			 default policies are to throw exceptions on arguments that
			 are out of domain or cause errors like numeric-overflow.
			 
			 Lacking try & catch blocks, the program will abort, whereas the
			 message below from the thrown exception will give some helpful
			 clues as to the cause of the problem.
			 */
			std::cout <<
			"\n""Message from thrown exception was:\n   " << e.what() << std::endl;
		}
		
	}
}

bool BundleFactory::spans_bad_intron(const ReadHit& read)
{

	const vector<CigarOp>& cig = read.cigar();
	
	size_t g_left = read.left();
	BadIntronTable::const_iterator itr = _bad_introns.find(read.ref_id());
	if (itr == _bad_introns.end())
		return false;
	
	const vector<AugmentedCuffOp>& bi = itr->second; 
	for (size_t i = 0; i < cig.size(); ++i)
	{
		assert(cig[i].length >= 0);
		switch(cig[i].opcode)
		{
			case MATCH:
				//ops.push_back(AugmentedCuffOp(CUFF_MATCH, g_left, cig[i].length));
				g_left += cig[i].length;
				break;
				
			case REF_SKIP:
			{	
				AugmentedCuffOp intron(CUFF_INTRON, g_left, cig[i].length);
				if (binary_search(bi.begin(), bi.end(), intron))
				{
					return true;
				}
				
				//ops.push_back(AugmentedCuffOp(CUFF_INTRON, g_left, cig[i].length));
				g_left += cig[i].length;
				break;
			}
				
			case SOFT_CLIP:
				g_left += cig[i].length;
				break;
                
            case HARD_CLIP:
				break;
            case INS:
                g_left -= cig[i].length;
                break;
            case DEL:
                g_left += cig[i].length;
                break;
			default:
				assert(false);
				break;
		}
	}
	
	return false;
}

void inspect_map(BundleFactory& bundle_factory,
                 BadIntronTable* bad_introns,
                 vector<LocusCount>& compatible_count_table,
                 vector<LocusCount>& total_count_table,
                 bool progress_bar,
                 bool show_stats)
{

	ProgressBar p_bar;
	if (progress_bar)
		p_bar = ProgressBar("Inspecting reads and determining fragment length distribution.",bundle_factory.ref_table().size());
	RefID last_chrom = 0;

	long double map_mass = 0.0;
    long double norm_map_mass = 0.0;
	
	int min_len = numeric_limits<int>::max();
	int max_len = def_max_frag_len;
	vector<double> frag_len_hist(def_max_frag_len+1,0);
	bool has_pairs = false;

	int num_bundles = 0;
	size_t total_hits = 0;
	size_t total_non_redundant_hits = 0;
	
	//To be used for quartile normalization
	vector<long double> mass_dist; 	
	
	// Store the maximum read length for "first" and "second" reads to report to user.
	int max_1 = 0;
	int max_2 = 0;
	
	shared_ptr<MultiReadTable> mrt(new MultiReadTable());
	
	while(true)
	{
		HitBundle* bundle_ptr = new HitBundle();
		
		bool valid_bundle = bundle_factory.next_bundle(*bundle_ptr);
		HitBundle& bundle = *bundle_ptr;

        if (emit_raw_counts && rc_stage)
           rc_write_counts(bundle_factory.ref_table().get_name(bundle.ref_id()), bundle);

        if (use_compat_mass) //only count hits that are compatible with ref transcripts
        {
            // Take raw mass even if bundle is "empty", since we could be out of refs
            // with remaining hits
            map_mass += bundle.compatible_mass();
//            if (lib_norm_method == QUARTILE && bundle.compatible_mass() > 0)
//            {
//                mass_dist.push_back(bundle.compatible_mass());
//            }
        }
        else if (use_total_mass) //use all raw mass
        { 
            
            // Take raw mass even if bundle is "empty", since we could be out of refs
            // with remaining hits
            map_mass += bundle.raw_mass();
//            if (lib_norm_method == QUARTILE && bundle.raw_mass() > 0)
//            {
//                mass_dist.push_back(bundle.raw_mass());
//            }
        }
        else
        {
            fprintf(stderr, "Error: hit counting scheme for normalization is not set!\n");
            assert(false);
            exit(1);
        }
		
		const RefSequenceTable& rt = bundle_factory.ref_table();
		const char* chrom = rt.get_name(bundle.ref_id());
		char bundle_label_buf[2048];
        if (chrom)
        {
            sprintf(bundle_label_buf, "%s:%d-%d", chrom, bundle.left(), bundle.right());
            verbose_msg("Inspecting bundle %s with %lu reads\n", bundle_label_buf, bundle.hits().size());
            
            vector<string> gene_ids;
            vector<string> gene_short_names;
            BOOST_FOREACH(shared_ptr<Scaffold> s, bundle.ref_scaffolds())
            {
                if (s->annotated_gene_id() != "")
                    gene_ids.push_back(s->annotated_gene_id());
                if (s->annotated_gene_name() != "")
                    gene_short_names.push_back(s->annotated_gene_name());
            }
            compatible_count_table.push_back(LocusCount(bundle_label_buf, floor(bundle.compatible_mass()), bundle.ref_scaffolds().size(), gene_ids, gene_short_names));
            total_count_table.push_back(LocusCount(bundle_label_buf, floor(bundle.raw_mass()), bundle.ref_scaffolds().size(), gene_ids, gene_short_names));
		}
        
        if (!valid_bundle)
		{
			delete bundle_ptr;
			break;
		}
		num_bundles++;
        
        if (progress_bar) 
        {
			double inc_amt = last_chrom == bundle.ref_id() ? 0.0 : 1.0;
			p_bar.update(bundle_label_buf, inc_amt);
			last_chrom = bundle.ref_id();
        }
        
        if (bad_introns != NULL)
		{
			identify_bad_splices(bundle, *bad_introns);
		}
		
		const vector<MateHit>& hits = bundle.non_redundant_hits();
		if (hits.empty())
		{
			delete bundle_ptr;
			continue;
		}
		
		list<pair<int, int> > open_ranges;
		int curr_range_start = hits[0].left();
		int curr_range_end = numeric_limits<int>::max();
		int next_range_start = -1;
		
		total_non_redundant_hits += bundle.non_redundant_hits().size();
		total_hits += bundle.hits().size();
		
		// This first loop calclates the map mass and finds ranges with no introns
		// Note that we are actually looking at non-redundant hits, which is why we use collapse_mass
		// This loop will also add multi-reads to the MultiReads table 
		for (size_t i = 0; i < hits.size(); ++i) 
		{
			assert(hits[i].left_alignment());
            
            // Add to table if multi-read
			if (hits[i].is_multi())
			{
				mrt->add_hit(hits[i]);
			}
			
			// Find left length
			int left_len = hits[i].left_alignment()->right()-hits[i].left_alignment()->left();
			min_len = min(min_len, left_len);
			if (!hits[i].left_alignment()->contains_splice())
            {
				if (hits[i].left_alignment()->is_first())
                    max_1 = max(max_1, left_len);
                else
                    max_2 = max(max_2, left_len);
            }
			
			// Find right length
			if (hits[i].right_alignment())
			{
				int right_len = hits[i].right_alignment()->right()-hits[i].right_alignment()->left();
				min_len = min(min_len, right_len);
				if (!hits[i].right_alignment()->contains_splice())
                {
                    if (hits[i].right_alignment()->is_first())
                        max_1 = max(max_1, right_len);
                    else
                        max_2 = max(max_2, right_len);
                }
                has_pairs = true;
			}
			
			// Find fragment length
			if (bundle.ref_scaffolds().size()==1 && hits[i].is_pair())
			// Annotation provided and single isoform gene
			{
				int start, end, mate_length;
				shared_ptr<Scaffold> scaff = bundle.ref_scaffolds()[0];
				if (scaff->map_frag(hits[i], start, end, mate_length))
				{
					if (mate_length >= min_len && mate_length <= max_len)
						frag_len_hist[mate_length] += hits[i].collapse_mass();
				}
			}
			else if (bundle.ref_scaffolds().empty())
			// No annotation provided.  Look for ranges.
			{
				if (hits[i].left() > curr_range_end)
				{
					if (curr_range_end - curr_range_start > max_len)
						open_ranges.push_back(make_pair(curr_range_start, curr_range_end));
					curr_range_start = next_range_start;
					curr_range_end = numeric_limits<int>::max();
				}
				if (hits[i].left_alignment()->contains_splice())
				{
					if (hits[i].left() - curr_range_start > max_len)
						open_ranges.push_back(make_pair(curr_range_start, hits[i].left()-1));
					curr_range_start = max(next_range_start, hits[i].left_alignment()->right());
				}
				if (hits[i].right_alignment() && hits[i].right_alignment()->contains_splice())
				{
					assert(hits[i].right_alignment()->left() >= hits[i].left());
					curr_range_end = min(curr_range_end, hits[i].right_alignment()->left()-1);
					next_range_start = max(next_range_start, hits[i].right());
				}
			}
		}
        
        if (bundle.ref_scaffolds().empty() && has_pairs) // No annotation provided
		{
			pair<int, int> curr_range(-1,-1);
			
			// This second loop uses the ranges found above to find the estimated frag length distribution
			// It also finds the minimum read length to use in the linear interpolation
			for (size_t i = 0; i < hits.size(); ++i)
			{
				if (hits[i].left() > curr_range.second && open_ranges.empty())
					break;
				
				if (hits[i].left() > curr_range.second)
				{
					curr_range = open_ranges.front();
					open_ranges.pop_front();
				}
				
				if (hits[i].left() >= curr_range.first && hits[i].right() <= curr_range.second && hits[i].is_pair())
				{
					int mate_len = hits[i].right()-hits[i].left();
					if (mate_len <= max_len)
						frag_len_hist[mate_len] += hits[i].collapse_mass();
				}
			}
		}
		
        open_ranges.clear();
		delete bundle_ptr;
	}
	
    norm_map_mass = map_mass;
    
//	if (lib_norm_method == QUARTILE && mass_dist.size() > 0)
//	{
//		sort(mass_dist.begin(),mass_dist.end());
//		int upper_quart_index = mass_dist.size() * 0.75;
//		norm_map_mass = mass_dist[upper_quart_index];
//	}

    if (bad_introns != NULL)
    {
        size_t alloced = 0;
        size_t used = 0;
        size_t num_introns = 0;
        for (BadIntronTable::const_iterator itr = bad_introns->begin();
             itr != bad_introns->end();
             ++itr)
        {
            alloced += itr->second.capacity() * sizeof(AugmentedCuffOp);
            used += itr->second.size() * sizeof(AugmentedCuffOp);
            num_introns += itr->second.size();
        }
        
        verbose_msg( "Bad intron table has %lu introns: (%lu alloc'd, %lu used)\n", num_introns, alloced, used);
    	verbose_msg( "Map has %lu hits, %lu are non-redundant\n", total_hits, total_non_redundant_hits);
    } 
    
	if (progress_bar)
		p_bar.complete();
	
	vector<double> frag_len_pdf(max_len+1, 0.0);
	vector<double> frag_len_cdf(max_len+1, 0.0);
    long double tot_count = accumulate(frag_len_hist.begin(), frag_len_hist.end(), 0.0 );
    bool empirical = false;
	
	if (user_provided_fld && has_pairs && tot_count >= 10000)
	{
		fprintf(stderr, "Warning: Overriding empirical fragment length distribution with user-specified parameters is not recommended.\n");
	}
	
	if (!has_pairs || tot_count < 10000)
	{
		if (has_pairs && !user_provided_fld)
		{
			fprintf(stderr, "Warning: Using default Gaussian distribution due to insufficient paired-end reads in open ranges.  It is recommended that correct parameters (--frag-len-mean and --frag-len-std-dev) be provided.\n");
		}
		tot_count = 0;
		normal frag_len_norm(def_frag_len_mean, def_frag_len_std_dev);
		max_len = def_frag_len_mean + 3*def_frag_len_std_dev;
		for(int i = min_len; i <= max_len; i++)
		{
			frag_len_hist[i] = cdf(frag_len_norm, i+0.5)-cdf(frag_len_norm, i-0.5);
			tot_count += frag_len_hist[i];
		}
	}
	else
	// Calculate the max frag length and interpolate all zeros between min read len and max frag len
	{	
		empirical = true;
		double curr_total = 0;
		size_t last_nonzero = min_len-1;
		for(size_t i = last_nonzero+1; i < frag_len_hist.size(); i++)
		{
			if (frag_len_hist[i] > 0)
			{
				if (last_nonzero != i-1)
				{
					double b = frag_len_hist[last_nonzero];
					double m = (frag_len_hist[i] - b)/(i-last_nonzero);
					for (size_t x = 1; x < i - last_nonzero; x++)
					{
						frag_len_hist[last_nonzero+x] = m * x + b;
						tot_count += frag_len_hist[last_nonzero+x];
						curr_total += frag_len_hist[last_nonzero+x];
					}	
				}
				last_nonzero = i;
			}
			
			curr_total += frag_len_hist[i];
			
			if (curr_total/tot_count > 0.9999)
			{
				max_len = i; 
				tot_count = curr_total;
				break;
			}
		}
	}
	
    double mean = 0.0;

    if (output_fld)
    {
        FILE* fhist = fopen(string(output_dir + "/frag_len_hist.csv").c_str(),"w");
        fprintf(fhist, "Length,Count\n");
        for(size_t i = 1; i < frag_len_hist.size(); i++)
        {
            fprintf(fhist, "%zu,%f\n", i, frag_len_hist[i]);
        }
        fclose(fhist);
    }

	// Convert histogram to pdf and cdf, calculate mean
	int frag_len_mode = 0;
	for(size_t i = min_len; i <= (size_t)max_len; i++)
	{
		frag_len_pdf[i] = frag_len_hist[i]/tot_count;
		frag_len_cdf[i] = frag_len_cdf[i-1] + frag_len_pdf[i];
        
		if (frag_len_pdf[i] > frag_len_pdf[frag_len_mode])
			frag_len_mode = i;
        mean += frag_len_pdf[i] * i;
	}
    
    double std_dev =  0.0;
    for(size_t i = 1; i < frag_len_hist.size(); i++)
    {
        std_dev += frag_len_pdf[i] * ((i - mean) * (i - mean));
    }
    
    std_dev = sqrt(std_dev);
	
	shared_ptr<ReadGroupProperties> rg_props = bundle_factory.read_group_properties();

    FLDSource source = DEFAULT;
    if (empirical)
    {
        source = LEARNED;
    }
    else if (user_provided_fld)
    {
        source = USER;
    }

	shared_ptr<EmpDist const> fld(new EmpDist(frag_len_pdf, frag_len_cdf, frag_len_mode, mean, std_dev, min_len, max_len, source));
	rg_props->multi_read_table(mrt);
	rg_props->frag_len_dist(fld);
	rg_props->normalized_map_mass(norm_map_mass);
    rg_props->total_map_mass(map_mass);

    if (show_stats)
    {
        fprintf(stderr, "> Map Properties:\n");
        //if (lib_norm_method == QUARTILE)
        //    fprintf(stderr, ">\tUpper Quartile: %.2Lf\n", norm_map_mass);
        fprintf(stderr, ">\tNormalized Map Mass: %.2Lf\n", norm_map_mass);
        fprintf(stderr, ">\tRaw Map Mass: %.2Lf\n", map_mass);
        if (corr_multi)
            fprintf(stderr,">\tNumber of Multi-Reads: %zu (with %zu total hits)\n", mrt->num_multireads(), mrt->num_multihits()); 
    //	if (has_pairs)
    //		fprintf(stderr, ">\tRead Type: %dbp x %dbp\n", max_1, max_2);
    //	else
    //		fprintf(stderr, ">\tRead Type: %dbp single-end\n", max(max_1,max_2));

        if (empirical)
        {
            fprintf(stderr, ">\tFragment Length Distribution: Empirical (learned)\n");
            fprintf(stderr, ">\t              Estimated Mean: %.2f\n", mean);
            fprintf(stderr, ">\t           Estimated Std Dev: %.2f\n", std_dev);
        }
        else
        {
            if (user_provided_fld)
            {
                fprintf(stderr, ">\tFragment Length Distribution: Truncated Gaussian (user-specified)\n");
            }
            else
            {
                fprintf(stderr, ">\tFragment Length Distribution: Truncated Gaussian (default)\n");
            }
            fprintf(stderr, ">\t              Default Mean: %d\n", def_frag_len_mean);
            fprintf(stderr, ">\t           Default Std Dev: %d\n", def_frag_len_std_dev);
        }
    }
	bundle_factory.num_bundles(num_bundles);
	bundle_factory.reset(); 
	return;
}
