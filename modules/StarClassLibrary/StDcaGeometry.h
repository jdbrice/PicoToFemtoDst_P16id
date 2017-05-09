/***************************************************************************
 *
 * 
 **************************************************************************/
#ifndef DCA_GEOMETRY_H
#define DCA_GEOMETRY_H
#include "TObject.h"
#include "StThreeVectorD.hh"
#include "StPhysicalHelixD.hh"
// #include "THelixTrack.h"

#include <iostream>
using namespace std;

class StDcaGeometry : public TObject {
public:
	StDcaGeometry();
	virtual ~StDcaGeometry();

	Int_t               charge()    const {return (mPti>0)? -1:1;}  // synchro with StiTrackNode charge definition
	Double_t            impact()    const {return mImp;}		  
	Double_t            curvature() const {return mCurv;}		  
	Double_t            psi()       const {return mPsi ;}		  
	Double_t            dipAngle()  const {return atan(mTan);}	  
	Double_t            tanDip()    const {return mTan ;}		  
	Double_t            pt()        const {return 1./fabs(mPti);}	  
	Double_t            z()         const {return mZ   ;}		  
	Double_t            hz()        const {return mCurv/mPti;}	  
	StThreeVectorD      origin()    const;         
	StThreeVectorD      momentum()  const;        
	StPhysicalHelixD    			helix()     const;
	// THelixTrack         thelix()    const;
	const float*        params()    const {return &mImp;}	  
	const float*        errMatrix() const {return &mImpImp;} 
	// void                GetXYZ(Double_t xyzp[6], Double_t CovXyzp[21]) const;
	virtual void        Print(Option_t *option = "") const;
	
	// Experts only set function
	void set(const Float_t pars[6], const Float_t errs[15]);
	void set(const Double_t pars[6], const Double_t errs[15]);

private:
	Char_t  mBeg[1];//!
	/// signed impact parameter; Signed in such a way that:
	///     x =  -impact*sin(Psi)
	///     y =   impact*cos(Psi)
	Float_t  mImp;
	///  Z-coordinate of this track (reference plane)
	Float_t  mZ;
	///  Psi angle of the track
	Float_t  mPsi;
	/// signed invert pt [sign = sign(-qB)]
	Float_t  mPti;
	/// tangent of the track momentum dip angle
	Float_t  mTan;
	/// signed curvature
	Float_t  mCurv;
	
	/// pars errors
	Float_t  mImpImp;
	Float_t  mZImp, mZZ;
	Float_t  mPsiImp, mPsiZ, mPsiPsi;
	Float_t  mPtiImp, mPtiZ, mPtiPsi, mPtiPti;
	Float_t  mTanImp, mTanZ, mTanPsi, mTanPti, mTanTan;
	Char_t   mEnd[1];//!
	
	// ClassDef(StDcaGeometry,1)
};
std::ostream&  operator<<(std::ostream& os, StDcaGeometry const & dca);

#endif
