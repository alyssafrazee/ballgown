#ifndef BUNDLES_H
#define BUNDLES_H
/*
 *  bundles.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 9/6/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <boost/bind.hpp>
#include <boost/random.hpp>
#include <vector>
#include <numeric>
#include "common.h"
#include "hits.h"
#include "scaffolds.h"
#include "gtf_tracking.h"
#include "progressbar.h"

extern bool rc_stage; //true lets HitBundle collect the raw counts

struct BundleStats
{		
	BundleStats() : 
	compatible(0), 
	uncollapsible(0),
	closure_edges(0),
	matched_edges(0)
	{}
	
	int compatible;
	int uncollapsible;
	int closure_edges;
	int matched_edges;
};


#define RC_MIN_EOVL 5
//Bundle raw count data
struct RC_Feature { //exon or intron of a reference transcript
	uint id; //feature id (>0)
	int l; int r; //genomic coordinates
	char strand;
	mutable uint rcount; //# reads covering this feature
	mutable uint ucount; //# uniquely mapped reads covering this feature
	mutable double mrcount; //multi-mapping-weighted counts
    //mutable vector<int> coverage; //per-base exon coverage data
	struct PCompare {
	 bool operator()(const RC_Feature* p1, const RC_Feature* p2) {
	 return (*p1 < *p2);
	 }
	};

	RC_Feature(int l0=0, int r0=0, char s='.', uint fid=0): id(fid), l(l0), r(r0),
		strand(s), rcount(0),ucount(0),mrcount(0) {
	if (l>r) { int t=l; l=r; r=t; }
	}

	bool operator<(const RC_Feature& o) const {
	 //if (id == o.id) return false;
	 if (l != o.l) return (l < o.l);
     if (r != o.r) return (r < o.r);
     if (strand == '.' || o.strand == '.') return false;
     if (strand != o.strand) return (strand < o.strand);
     return false;
	 }
	bool operator==(const RC_Feature& o) const {
	 //if (id == o.id) return true;
	 return (l==o.l && r==o.r &&
		 (strand == o.strand || strand == '.' || o.strand == '.'));
	 }
	bool strand_compatible(const RC_Feature& o) const {
		 return (strand == '.' || o.strand == '.' || strand == o.strand);
	}
	//WARNING: the overlap checks IGNORE strand!
	bool overlap(int hl, int hr) const {
	  if (hl>hr) { int t=hl; hl=hr; hr=t; }
      return (l<=hr && r<=hl);
	  }
	bool overlap(int hl, int hr, int minovl) const {
	  if (hl>hr) { int t=hl; hl=hr; hr=t; }
      hl+=minovl;hr-=minovl;
      return (l<=hr && r<=hl);
	  }
	uint ovlen(int hl, int hr) const {
     if (hl>hr) { int t=hl; hl=hr; hr=t; }
     if (l<hl) {
        if (hl>r) return 0;
        return (hr>r) ? r-hl+1 : hr-hl+1;
        }
       else { //hl<=l
        if (l>hr) return 0;
        return (hr<r)? hr-l+1 : r-l+1;
        }
	 }
};


typedef set<const RC_Feature*, RC_Feature::PCompare> RC_FeatPtrSet;
typedef set<RC_Feature>::iterator RC_FeatIt;
typedef map<uint, set<uint> > RC_Map2Set;
typedef map<uint, set<uint> >::iterator RC_Map2SetIt;
struct RC_Seg { //just a genomic interval holder
	int l;
	int r;
	RC_Seg(int l0=0, int r0=0):l(l0), r(r0) { }
};

struct RC_ScaffData {
	const Scaffold* scaff;
	string t_name; //original GFF ID for the transcript
	int l;
	int r;
	char strand;
	int t_id;
	int num_exons;
	int eff_len;
	mutable double cov;
	mutable double fpkm;
	//other mutable fields here, to be updated by rc_update_scaff()
	RC_ScaffData(const Scaffold* s=NULL):scaff(s), t_name(), l(0), r(0),
		strand('.'), t_id(0), num_exons(0), eff_len(0), cov(0), fpkm(0) {
	  if (scaff==NULL || scaff->rc_id_data()==NULL) return;
	  RC_ScaffIds& sdata = *(scaff->rc_id_data());
	  t_id = sdata.t_id;
	  t_name=scaff->annotated_trans_id();
	  strand=sdata.strand;
	  l=scaff->left();
	  r=scaff->right();
	  num_exons=sdata.exons.size();
	  for (size_t i=0;i<sdata.exons.size();++i) {
		RC_ScaffSeg& exon = sdata.exons[i];
		eff_len+=exon.r-exon.l;
	  }
	}

    bool operator<(const RC_ScaffData& o) const {
    	if (l != o.l) return (l < o.l);
    	if (r != o.r) return (r < o.r);
    	if (strand != o.strand) return (strand < o.strand);
	    return (t_name < o.t_name);
		return false;
    }
};

FILE* rc_fwopen(const char* fname);
int rc_cov_inc(int i);
class RC_MultiCovInc {
	float fcov;
  public:
	RC_MultiCovInc(int numhits):fcov(1.0) {
	 if (numhits>1) fcov=1/(float)numhits;
	}
	float operator()(const float& v) {
	  return (v+fcov);
	}
};
struct RC_BundleData {
 int lmin;
 int rmax;
 set<RC_ScaffData> tdata;
 map<uint, set<uint> > e2t; //mapping exon ID to transcript IDs
 map<uint, set<uint> > i2t; //mapping intron ID to transcript IDs
 set<RC_Feature> exons; //exons by their start coordinate
 set<RC_Feature> introns; //introns by their start coordinate
 RC_FeatIt xcache; //cache the first exon overlapping xcache_pos to speed up exon-overlap queries (findExons())
 int xcache_pos; // left coordinate of last cached exon overlap query (findExons())
 // -- output files
 FILE* ftdata; //t_data
 FILE* fedata; //e_data
 FILE* fidata; //i_data
 FILE* fe2t;   //e2t
 FILE* fi2t;   //i2t
 vector<float> f_mcov; //coverage data, multi-map aware, per strand
 vector<int> f_cov;
 vector<float> r_mcov; //coverage data on the reverse strand
 vector<int> r_cov;
 //
 RC_BundleData(int bundle_l=0, int bundle_r=0):lmin(bundle_l),rmax(bundle_r),
	 tdata(), e2t(), i2t(), exons(), introns(), xcache(exons.end()),
	 xcache_pos(0), ftdata(NULL), fedata(NULL), fidata(NULL),
	 fe2t(NULL), fi2t(NULL) { }
 void setupFiles(FILE* &f_tdata, FILE* &f_edata, FILE* &f_idata,
	            FILE* &f_e2t, FILE* &f_i2t) {
   if (f_tdata == NULL) {
 	//first call, create the files
	 f_tdata = rc_fwopen("t_data");
	 fprintf(f_tdata, "t_id\tchr\tstrand\tstart\tend\tt_name\tnum_exons\tlength\tgene_id\tgene_name\tcov\tFPKM\n");
	 f_edata = rc_fwopen("e_data");
     fprintf(f_edata, "e_id\tchr\tstrand\tstart\tend\trcount\tucount\tmrcount\tcov\tcov_sd\tmcov\tmcov_sd\n");
	 f_idata = rc_fwopen("i_data");
     fprintf(f_idata, "i_id\tchr\tstrand\tstart\tend\trcount\tucount\tmrcount\n");
     f_e2t = rc_fwopen("e2t");
     fprintf(f_e2t,  "e_id\tt_id\n");
     f_i2t = rc_fwopen("i2t");
     fprintf(f_i2t,  "i_id\tt_id\n");
   }
  ftdata=f_tdata;
  fedata=f_edata;
  fidata=f_idata;
  fe2t=f_e2t;
  fi2t=f_i2t;
 }
 void addFeature(uint t_id, int l, int r, char strand, uint f_id, set<RC_Feature>& fset,
	                         map<uint, set<uint> >& f2t) {
   RC_Feature feat(l, r, strand, f_id);
   pair<RC_FeatIt, bool> in = fset.insert(feat);
   //if (!in.second) { //existing f_id
   // f_id=in.first->id;
   //}
   set<uint> tset;
   tset.insert(t_id);
   pair<RC_Map2SetIt, bool> mapin=f2t.insert(pair<uint, set<uint> >(f_id, tset));
   if (!mapin.second) {
	 //existing f_id
	 (*mapin.first).second.insert(t_id);
   }
  }

 void addScaffold(Scaffold& ps) {
   if (!ps.rc_id_data()) return;
   RC_ScaffIds& sdata = *(ps.rc_id_data());
   RC_ScaffData scaffdata(&ps);
   tdata.insert(scaffdata);
   if (lmin==0 || lmin>ps.left()) lmin=ps.left();
   if (rmax==0 || rmax<ps.right()) rmax=ps.right();
   for (vector<RC_ScaffSeg>::iterator it=sdata.exons.begin();it!=sdata.exons.end();++it) {
	 addFeature(sdata.t_id, it->l, it->r, sdata.strand, it->id, exons, e2t);
   }
   //store introns too
   for (vector<RC_ScaffSeg>::iterator it=sdata.introns.begin();it!=sdata.introns.end();++it) {
	 addFeature(sdata.t_id, it->l, it->r, sdata.strand, it->id, introns, i2t);
   }
 }
 void setupCov() {
   //to be called after all reference transcripts were added
   assert(rmax>lmin);
   int blen=rmax-lmin+1;
   f_cov.resize(blen, 0);
   r_cov.resize(blen, 0);
   f_mcov.resize(blen, 0.0);
   r_mcov.resize(blen, 0.0);
 }

 void updateCov(char strand, int numhits, int gpos, int glen) {
  if (gpos>rmax || gpos+glen<lmin) return; //no overlap with bundle
  if (gpos<lmin) { //read overlap begins before the bundle start (!)
	int gadj=lmin-gpos;
	gpos+=gadj;
	glen-=gadj;
  }
  if (gpos+glen>rmax) {
	glen=rmax-gpos;
  }
  if (glen<=0) return; //no overlap (shouldn't get here)
  int goffs=gpos-lmin;
  if (goffs<0) return; //should NOT be here!
  if (strand=='.' || strand=='+') {
	transform(f_cov.begin()+goffs, f_cov.begin()+goffs+glen,
		        f_cov.begin()+goffs, rc_cov_inc);
    transform(f_mcov.begin()+goffs, f_mcov.begin()+goffs+glen,
    	        f_mcov.begin()+goffs, RC_MultiCovInc(numhits));
  }
  if (strand=='.' || strand=='-') {
	transform(r_cov.begin()+goffs, r_cov.begin()+goffs+glen,
		        r_cov.begin()+goffs, rc_cov_inc);
    transform(r_mcov.begin()+goffs, r_mcov.begin()+goffs+glen,
    	        r_mcov.begin()+goffs, RC_MultiCovInc(numhits));
  }

 }

 RC_FeatPtrSet findExons(int hl, int hr, char strand='.', bool update_cache=true) {
  //returns exons overlapping given interval hl-hr
   RC_FeatPtrSet ovlex; //return set
   RC_Feature q(hl, hr);
   RC_FeatIt xstart=exons.begin();
   bool no_cache=(xcache_pos==0 || xcache_pos>hl);
   if (no_cache) {
	   if (update_cache) {
	      xcache=exons.end();
	      xcache_pos=0;
	   }
   }
   else xstart=xcache; //must have a valid value
   bool upd_cache(update_cache);
   RC_FeatIt last_checked_exon(exons.end());
   for (RC_FeatIt p=xstart;p != exons.end();++p) {
     last_checked_exon=p;
	 if (p->l > hr) break;
	 if (hl > p->r) continue;
     //exon overlap
     if (upd_cache) {
         //cache first overlap
		 xcache=p;
		 upd_cache=false;
     }
     if (strand!='.' && strand!=p->strand) continue;
	 ovlex.insert(&(*p));
   }
   if (update_cache) {
	 if (upd_cache) xcache=last_checked_exon; //there was no overlap found
	 xcache_pos=hl;
     }
   return ovlex;
  }

 RC_FeatIt findIntron(int hl, int hr, char strand) {
   RC_FeatIt ri=introns.find(RC_Feature(hl, hr, strand));
   return ri;
 }
}; //struct BundleRC_Data


typedef map<RefID, vector<AugmentedCuffOp> > BadIntronTable;


/*******************************************************************************
 HitBundle is a set of MateHit objects that, were you to look at the interval
 graph of their spanning intervals in genomic coordinates, you'd see a single
 connected component. Note that bundles do not correspond to single transcripts,
 or even single genes
 *******************************************************************************/
class HitBundle
{
private:
    HitBundle(const HitBundle& rhs) {} 
public:
	HitBundle() 
    : _leftmost(INT_MAX), _rightmost(-1), _final(false), _id(++_next_id), _ref_id(0), _raw_mass(0.0), _num_replicates(1), 
    _compatible_mass(0.0),rc_data(NULL) {}
	
    ~HitBundle()
    {
        delete rc_data;
        vector<shared_ptr<Scaffold> >& bundle_ref_scaffs = ref_scaffolds();
        BOOST_FOREACH(shared_ptr<Scaffold>& ref_scaff, bundle_ref_scaffs)
        {
            // This bundle and the factory that actually owns the ref_mRNAs
            // are the only objects that should have access to these scaffolds
            // so if the use count is 2, we can clear these guys.
			// Updated to 3 since the bias learner now uses them.
            if (ref_scaff.use_count() <= 3)
            {
                ref_scaff->clear_hits();
            }
            else if (ref_scaff->mate_hits().size() > 0)
            {
                fprintf(stderr, "Warning: bundle %d-%d shared reference scaffolds with others.  Possible soft memory leak.\n", left(), right());
            }
        }   
        
        BOOST_FOREACH (MateHit& hit, _hits)
		{
			delete hit.left_alignment();
			delete hit.right_alignment();
		}
        
        for(OpenMates::iterator itr = _open_mates.begin(); itr != _open_mates.end(); ++itr)
		{
			BOOST_FOREACH (MateHit& hit,  itr->second)
            {
                delete hit.left_alignment();
                delete hit.right_alignment();
            }
		}
		
    }
	int left()   const { return _leftmost;  }
	int right()  const { return _rightmost; }
	int length() const { return _rightmost - _leftmost; }
	
	// Returns true if the hit was added successfully.
	bool add_hit(const MateHit& hit);
    
	void rc_init() {
	  if (rc_data==NULL) {
		rc_data = new RC_BundleData(this->_leftmost, this->_rightmost);
	    }
	  }
	//after reference annotation was loaded:
	void rc_finalize_refs() {
      if (rc_data==NULL) return;
      rc_data->setupCov();
	}
	//store ref transcripts' features into rc_store
    void rc_store_t(Scaffold& scaff);
	//check read mapping against features of ref transcripts:
    void rc_count_hit(const ReadHit* bh);

	// This is to keep track of mass of all hits, including
	// thosethat are not added to any bundle
	// but are skipped during the creation of this bundle
	void add_raw_mass(double rm) { _raw_mass += rm; }
	void rem_raw_mass(double rm) { _raw_mass -= rm; }
	double raw_mass() { return _raw_mass; }
    
    double compatible_mass() const 
    {
        return _compatible_mass;
    }
    
	void clear_hits() 
    {
        _hits.clear(); 
        _non_redundant.clear();
        vector<shared_ptr<Scaffold> >& bundle_ref_scaffs = ref_scaffolds();
        BOOST_FOREACH(shared_ptr<Scaffold>& ref_scaff, bundle_ref_scaffs)
        {
            if (ref_scaff.use_count() <= 3)
            {
                ref_scaff->clear_hits();
            }
            else 
            {
                fprintf(stderr, "Warning: bundle %d-%d shared reference scaffolds with others.  Possible soft memory leak.\n", left(), right());
            }
        } 
    }
	
    const std::vector<MateHit>& hits() const { return _hits; } 
	const std::vector<MateHit>& non_redundant_hits() const { return _non_redundant; } 
	
	RefID ref_id()  const {return _ref_id; }
	
    RC_BundleData* rcdata() const { return rc_data; }

	int id() const { return _id; }
	
	void add_ref_scaffold(shared_ptr<Scaffold> scaff)
	{
		if (scaff->left() < _leftmost)
			_leftmost = scaff->left();
		if (scaff->right() > _rightmost)
			_rightmost = scaff->right();
		_ref_scaffs.push_back(scaff);
        _ref_id = scaff->ref_id();
        //emit_raw_counts: store exons and introns for easy hit checking
    	if (emit_raw_counts && rc_stage) {
          rc_store_t(*scaff);
    	}

	}
	
	vector<shared_ptr<Scaffold> >& ref_scaffolds() { return _ref_scaffs; }
	
	// Adds a Bowtie hit to the open hits buffer.  The Bundle will handle turning
	// the Bowtie hit into a properly mated Cufflinks hit record
	bool add_open_hit(shared_ptr<ReadGroupProperties const> rg_props,
                      const ReadHit* bh,
					  bool expand_by = true);
	
	// Commits any mates still open as singleton hits
	void finalize_open_mates();
	
	// Sorts the hits, and performs other functions needed to prepare this 
	// bundle for assembly and quantitation
	void finalize(bool is_combined=false);
	
	void remove_hitless_scaffolds();

	void collapse_hits();
	
    int num_replicates() const { return _num_replicates; }
    
	double mass() const
	{
		double mass = 0;
		for(size_t i = 0; i < _non_redundant.size(); i++)
		{
			mass += _non_redundant[i].collapse_mass();
		}
		return mass;
	}
	
    static void combine(const vector<HitBundle*>& in_bundles,
                        HitBundle& out_bundle);
    
private:
    int _leftmost;
	int _rightmost;
	std::vector<MateHit> _hits;
	std::vector<MateHit> _non_redundant;
	std::vector<shared_ptr<Scaffold> > _ref_scaffs; // user-supplied reference annotations overlapping the bundle
	bool _final;
	int _id;
    RefID _ref_id;
	double _raw_mass;

	
	static int _next_id;
	
	typedef map<int, list<MateHit> > OpenMates;
	OpenMates _open_mates;
    int _num_replicates;
    double _compatible_mass;

	RC_BundleData* rc_data;
};

void load_ref_rnas(FILE* ref_mRNA_file, 
				   RefSequenceTable& rt,
				   vector<shared_ptr<Scaffold> >& ref_mRNAs,
				   bool loadSeqs=false,
				   bool loadFPKM=false);

class BundleFactory
{
public:
    
	BundleFactory(shared_ptr<HitFactory> fac, BundleMode bm)
	: _hit_fac(fac), _bundle_mode(bm), _prev_pos(0), _prev_ref_id(0), _curr_bundle(0),  _zeroone(rng)
	{
		_rg_props = shared_ptr<ReadGroupProperties>(new ReadGroupProperties(fac->read_group_properties()));
        
        
       
	}

    bool bundles_remain()  
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        return _curr_bundle < num_bundles();
    }
    
	bool next_bundle(HitBundle& bundle_out);
	bool next_bundle_hit_driven(HitBundle& bundle_out);
	bool next_bundle_ref_driven(HitBundle& bundle_out);
	bool next_bundle_ref_guided(HitBundle& bundle_out);

    
    RefSequenceTable& ref_table() { return _hit_fac->ref_table(); }
    
    // Not available until after inspect_bundle
    int num_bundles() const { return _num_bundles; }
    void num_bundles(int n) { _num_bundles = n; }
    
	void reset() 
	{ 
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        _curr_bundle = 0;
		//rewind(hit_file); 
		_hit_fac->reset();
		next_ref_scaff = ref_mRNAs.begin(); 
        next_mask_scaff = mask_gtf_recs.begin();
        
        BOOST_FOREACH(shared_ptr<Scaffold> ref_scaff, ref_mRNAs)
        {
            ref_scaff->clear_hits();
        }
        
        _prev_pos = 0;
        _prev_ref_id = 0;
	}
	
    // This function NEEDS to deep copy the ref_mRNAs, otherwise cuffdiff'd
    // samples will clobber each other
    void set_ref_rnas(const vector<shared_ptr<Scaffold> >& mRNAs)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        ref_mRNAs.clear();
        for (vector<shared_ptr<Scaffold> >::const_iterator i = mRNAs.begin(); i < mRNAs.end(); ++i)
        {
            ref_mRNAs.push_back(shared_ptr<Scaffold>(new Scaffold(**i)));
        }
        
        RefID last_id = 0;
        for (vector<shared_ptr<Scaffold> >::iterator i = ref_mRNAs.begin(); i < ref_mRNAs.end(); ++i)
        {
            if ((*i)->ref_id() != last_id)
            {
                _ref_scaff_offsets.push_back(make_pair((*i)->ref_id(), i));
            }
            last_id = (*i)->ref_id();
        }
        
        next_ref_scaff = ref_mRNAs.begin();
    }
    
    void set_mask_rnas(const vector<shared_ptr<Scaffold> >& masks)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        mask_gtf_recs = masks;
        RefID last_id = 0;
        for (vector<shared_ptr<Scaffold> >::iterator i = mask_gtf_recs.begin(); i < mask_gtf_recs.end(); ++i)
        {
            if ((*i)->ref_id() != last_id)
            {
                _mask_scaff_offsets.push_back(make_pair((*i)->ref_id(), i));
            }
            last_id = (*i)->ref_id();
        }
        
        next_mask_scaff = mask_gtf_recs.begin();
    }
        
	void bad_intron_table(const BadIntronTable& bad_introns) 
	{ 
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
		_bad_introns = bad_introns;
	}
    
    void read_group_properties(shared_ptr<ReadGroupProperties> rg)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        _rg_props = rg;
    }
    
    shared_ptr<ReadGroupProperties> read_group_properties()
    {
        return _rg_props;
    }
	
	bool spans_bad_intron(const ReadHit& read);
	
private:
	
	bool _expand_by_hits(HitBundle& bundle);
	bool _expand_by_refs(HitBundle& bundle);
	
	shared_ptr<HitFactory> _hit_fac;
    
	vector<shared_ptr<Scaffold> > ref_mRNAs;
	//FILE* ref_mRNA_file;
	vector<pair<RefID, vector<shared_ptr<Scaffold> >::iterator> > _ref_scaff_offsets;
	vector<shared_ptr<Scaffold> >::iterator next_ref_scaff;
    
    vector<shared_ptr<Scaffold> > mask_gtf_recs;
	//FILE* mask_file;
	vector<pair<RefID, vector<shared_ptr<Scaffold> >::iterator> > _mask_scaff_offsets;
	vector<shared_ptr<Scaffold> >::iterator next_mask_scaff;
	
	BadIntronTable _bad_introns;
    
    shared_ptr<ReadGroupProperties> _rg_props;
    
	// Sets nva to point to the next valid alignment
	// Returns the mass of any alignments that are seen, valid or not
    double next_valid_alignment(const ReadHit*& nva);
	
	// Backs up the factory to before the last valid alignment
	// and returns the mass of that alignment (rh)
	double rewind_hit(const ReadHit* rh);

    BundleMode _bundle_mode;
    int _prev_pos;
    RefID _prev_ref_id;
    int _num_bundles;
    int _curr_bundle;
    
    boost::mt19937 rng;
    boost::uniform_01<boost::mt19937> _zeroone;
    
#if ENABLE_THREADS    
    boost::mutex _factory_lock;
#endif
};

//raw counts
void rc_update_tdata(HitBundle& bundle, const Scaffold& scaff, double cov, double fpkm);

void rc_write_counts(const char* refname, HitBundle& bundle);


void identify_bad_splices(const HitBundle& bundle, 
						  BadIntronTable& bad_splice_ops);

void inspect_map(BundleFactory& bundle_factory,
                 BadIntronTable* bad_introns,
                 vector<LocusCount>& compatible_count_table,
                 vector<LocusCount>& total_count_table,
                 bool progress_bar = true,
                 bool show_stats = true);

#endif
