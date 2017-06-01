#ifndef FEMTO_DST_WRITER_H
#define FEMTO_DST_WRITER_H

#include "PicoDstSkimmer/PicoDstSkimmer.h"

#include "ProductionUtils/RunMapFactory.h"

#include "PicoDstP16id/StPicoMtdPidTraits.h"

#include "FemtoDstFormat/BranchWriter.h"
#include "FemtoDstFormat/TClonesArrayWriter.h"
#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrack.h"
#include "FemtoDstFormat/FemtoMtdPidTraits.h"
#include "FemtoDstFormat/FemtoBTofPidTraits.h"
#include "FemtoDstFormat/FemtoTrackHelix.h"

#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"

#include "TTree.h"
#include "TAxis.h"

class FemtoDstWriter : public PicoDstSkimmer
{
protected:
	FemtoEvent _event;
	RunMapFactory rmf;

	TTree * tree;
	BranchWriter<FemtoEvent> _wEvent;
	TClonesArrayWriter<FemtoTrack> _wTrack;
	TClonesArrayWriter<FemtoMtdPidTraits> _wMtdPid;
	TClonesArrayWriter<FemtoBTofPidTraits> _wBTofPid;
	TClonesArrayWriter<FemtoTrackHelix>    _wHelix;
	FemtoTrack _track;
	FemtoMtdPidTraits _mtdPid;
	FemtoBTofPidTraits _btofPid;
	FemtoTrackHelix _helix;

	StRefMultCorr *rmc = nullptr;

	bool require_mtdPid = false;
	bool require_btofPid = false;

	double max_pT = 1000;

public:
	virtual const char* classname() const { return "FemtoDstWriter"; }
	FemtoDstWriter() {}
	~FemtoDstWriter() {}

	virtual void initialize(){
		PicoDstSkimmer::initialize();

		book->cd();
		tree = new TTree( "FemtoDst", "FemtoDst" );
		_wEvent.createBranch( tree, "Event" );
		_wTrack.createBranch( tree, "Tracks" );
		// _wMtdPid.createBranch( tree, "MtdPidTraits" );
		_wBTofPid.createBranch( tree, "BTofPidTraits" );
		_wHelix.createBranch( tree, "Helices" );


		if ( "P16id" == config[nodePath + ".Require:rmc"] ){
			LOG_F( INFO, "Using RefMultCorr for P16id" );
			rmc = CentralityMaker::instance()->getgRefMultCorr_P16id();
			rmc->init(15075008);
			rmc->setVzForWeight(6, -6.0, 6.0);
			rmc->readScaleForWeight("config/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");
		} else {
			LOG_F( INFO, "Using RefMultCorr for MTD (P16id)" );
			rmc = CentralityMaker::instance()->getgRefMultCorr_VpdMBnoVtx();
		}
		

		require_mtdPid = config.getBool( nodePath + ".Require:mtd", false );
		require_btofPid = config.getBool( nodePath + ".Require:btof", false );
		max_pT = config.getDouble( nodePath + ".Require:max_pT", max_pT );

	}
protected:	

	virtual void preEventLoop(){
		TreeAnalyzer::preEventLoop();

	}
	virtual bool keepEvent( StPicoEvent *event ){

		book->fill( "events", "All" );
		// if ( false == event->isDiMuon() )
		// 	return false;

		book->fill( "events", "Trigger" );

		rmc->init( event->runId() );
		rmc->initEvent( 
			event->grefMult(), 
			event->primaryVertex().z(), 
			event->ZDCx()
		);

		auto vtx = event->primaryVertex();
		float deltaVz = event->vzVpd() - vtx.z();

		book->fill( "vtx_z", vtx.z() );
		book->fill( "delta_z", deltaVz );

		if ( fabs( vtx.z() ) > 100 )
			return false;
		book->fill( "events", "vtx" );

		if ( fabs(deltaVz) > 3.0 )
			return false;
		book->fill( "events", "vtx_delta" );


		if ( rmc->getCentralityBin16() < 0 || rmc->getCentralityBin16() > 16 )
			return false;
		book->fill( "events", "cent80" );


		book->fill( "pass_vtx_z", vtx.z() );
		book->fill( "pass_delta_z", deltaVz );

		return true;
	}

	virtual void analyzeEvent(){
		StPicoEvent *event = _rEvent.get( 0 );
		if ( nullptr == event ){
			LOG_F( ERROR, "NULL EVENT" );
		}

		if ( false == keepEvent( event ) )
			return;

		auto vtx = event->primaryVertex();

		_event.mPrimaryVertex_mX1 = vtx.x();
		_event.mPrimaryVertex_mX2 = vtx.y();
		_event.mPrimaryVertex_mX3 = vtx.z();

		_event.mRunId    = event->runId();
		_event.mRunIndex = rmf.indexForRun( _event.mRunId );
		_event.mEventId  = event->eventId();
		// _event.mTriggerWordMtd = event->triggerWordMtd();
		_event.mGRefMult = event->grefMult();
		_event.mBin16 = rmc->getCentralityBin16();
		_event.mWeight = rmc->getWeight();

		_wEvent.set( _event );

		_wTrack.reset();
		// _wMtdPid.reset();
		_wBTofPid.reset();
		_wHelix.reset();
		size_t nTracks = _rTrack.N();
		size_t nMtdTracks = 0;
		size_t nBTofTracks = 0;
		for ( size_t i = 0; i < nTracks; i++ ){
			StPicoTrack * track = _rTrack.get( i );

			fillTrack( nMtdTracks, track, event );
			if ( fabs(_track.mPt) > 0.01 && fabs(_track.mPt) < max_pT ){
				

				if ( require_mtdPid  && _track.mMtdPidTraitsIndex < 0 ) continue;
				if ( require_btofPid && _track.mBTofPidTraitsIndex < 0 ) {continue;}

				_wTrack.add( _track );
				_wHelix.add( _helix );
				// if ( _track.mMtdPidTraitsIndex >= 0 && require_mtdPid){
				// 	_wMtdPid.add( _mtdPid );
				// 	nMtdTracks++;
				// }

				if ( _track.mBTofPidTraitsIndex >= 0 ){
					nBTofTracks++;
					_wBTofPid.add( _btofPid );
				}
			}// mPt > 0.01
		} // loop on tracks

		if ( nMtdTracks >= 1 )
			book->fill( "events", "gte1_mtd" );
		if ( nMtdTracks >= 2 )
			book->fill( "events", "gte2_mtd" );

		tree->Fill();
	}

	virtual void fillTrack( size_t i, StPicoTrack *track, StPicoEvent *event ){
		// LOG_F( INFO, "Filling %lu", i );
		_track.reset();

		_track.mId  = track->id();
		auto pMom   = track->pMom();
		_track.mPt  = pMom.perp();
		_track.mEta = pMom.pseudoRapidity();
		_track.mPhi = pMom.phi();
		
		_track.mNHitsFit  = track->nHitsFit() * track->charge();
		_track.mNHitsMax  = track->nHitsMax();
		_track.mNHitsDedx = track->nHitsDedx();

		_track.nSigmaElectron ( track->nSigmaElectron ( ) );
		_track.nSigmaPion     ( track->nSigmaPion     ( ) );
		_track.nSigmaKaon     ( track->nSigmaKaon     ( ) );
		_track.nSigmaProton   ( track->nSigmaProton   ( ) );

		double dca_full_prec = track->globalDCA( event->bField(), event->primaryVertex() );
		_track.gDCA( dca_full_prec );

		_track.dEdx( track->dEdx() );


		fillTrackHelix( i, track, dca_full_prec );

		if ( track->mtdPidTraitsIndex() >= 0 ){
			fillMtdPid( i, track );
		}

		if ( track->bTofPidTraitsIndex() >= 0 ){
			fillBTofPid( i, track );
		}
		
	}


	virtual void fillMtdPid( size_t i, StPicoTrack *track ){
		_mtdPid.reset();

		size_t index = _wMtdPid.N();
		StPicoMtdPidTraits *mtdPid = _rMtdPid.get( track->mtdPidTraitsIndex() );
		_mtdPid.mMatchFlag = mtdPid->matchFlag();
		_mtdPid.mDeltaY = mtdPid->deltaY();
		_mtdPid.mDeltaZ = mtdPid->deltaZ();
		_mtdPid.mDeltaTimeOfFlight = mtdPid->deltaTimeOfFlight();
		_mtdPid.mMtdHitChan = mtdPid->mMtdHitChan;


		size_t nMtdHits = _rMtdHit.N();
		for ( size_t iMtdHit = 0; iMtdHit < nMtdHits; iMtdHit++ ){
			StPicoMtdHit *mtdHit = _rMtdHit.get( iMtdHit );
			
			if ( nullptr != mtdHit && mtdHit->gChannel() == _mtdPid.mMtdHitChan ){
				_mtdPid.mTriggerFlag = mtdHit->triggerFlag();
				break;
			}
		} // loop iMtdHit

		_track.mMtdPidTraitsIndex = index;

	}


	virtual void fillBTofPid( size_t i, StPicoTrack *track ){
		_btofPid.reset();
		StPicoBTofPidTraits *btofPid = _rBTofPid.get( track->bTofPidTraitsIndex() );

		_btofPid.beta( btofPid->btofBeta() );
		_btofPid.matchFlag( btofPid->btofMatchFlag() );
		_btofPid.yLocal( btofPid->btofYLocal() );
		_btofPid.zLocal( btofPid->btofZLocal() );

		size_t index = _wBTofPid.N();
		_track.mBTofPidTraitsIndex = index;

	}

	virtual void fillTrackHelix( size_t i, StPicoTrack *track, double dca ){
		_helix.reset();

		_helix.mPar[0] = track->params()[0];
		_helix.mPar[1] = track->params()[1];
		_helix.mPar[2] = track->params()[2];
		_helix.mPar[3] = track->params()[3];
		_helix.mPar[4] = track->params()[4];
		_helix.mPar[5] = track->params()[5];
		_helix.mMap0 = track->map0();
		_helix.mMap1 = track->map1();
		_helix.mDCA = dca;

		size_t index = _wHelix.N();
		_track.mHelixIndex = index;
	}

};


#endif