/***************************************************************************
 *
 *
 *
 *
 *
 **************************************************************************/
#include "StDcaGeometry.h"
#include "TCernLib.h"
#include "TMatrix.h"
#include "TMath.h"

// ClassImp(StDcaGeometry);
    
static const char rcsid[] = "$Id: StDcaGeometry.cxx,v 2.11 2013/11/13 21:35:48 fisyak Exp $";

// ClassImp(StDcaGeometry);

StDcaGeometry::StDcaGeometry()
{
    memset(mBeg,0,mEnd-mBeg+1);
}

StDcaGeometry::~StDcaGeometry() {/* noop */}

StThreeVectorD StDcaGeometry::origin() const
{
    float x = -mImp*sin(mPsi);
    float y =  mImp*cos(mPsi);
    return StThreeVectorD(x,y,mZ);
}

StThreeVectorD StDcaGeometry::momentum() const
{
    float ptt = pt();
    float x   = ptt*cos(mPsi);
    float y   = ptt*sin(mPsi);
    float z   = ptt*mTan;
    return StThreeVectorD(x,y,z);
}

void StDcaGeometry::set(const float pars[7],const float errs[15])
{
    if (pars) memcpy(&mImp   ,pars,sizeof(float)*6 );
    if (errs) memcpy(&mImpImp,errs,sizeof(float)*15);
}
void StDcaGeometry::set(const double pars[7],const double errs[15])
{
  if (pars) TCL::ucopy(pars, &mImp, 6);
  if (errs) TCL::ucopy(errs, &mImpImp, 15);
}

StPhysicalHelixD StDcaGeometry::helix() const
{
    //    double curvature = fabs(mCurv);
    int  h = (mCurv>=0) ? 1:-1;
    
    double phase = mPsi-h*M_PI/2;
    
    return StPhysicalHelixD(fabs(mCurv),   // 1/cm
                    atan(mTan),    // radian
                    phase,         // radian
                    origin(),      // cm
                    h);
}


ostream&  operator<<(ostream& os, const StDcaGeometry& dca) {
  const Float_t *errMx =  dca.errMatrix();
  return os << Form("Dca: imp %7.2f +/-%7.2f,Z:%7.2f +/-%7.2f,psi:%7.2f +/-%7.2f,pT/q:%7.2f +/-%6.1f%%,TanL:%8.3f +/-%8.3f",
		    dca.impact(),    (errMx[0] >= 0)  ? TMath::Sqrt(errMx[0]) : -13,
		    dca.z(),         (errMx[2] >= 0)  ? TMath::Sqrt(errMx[2]) : -13,
		    dca.psi(),       (errMx[5] >= 0)  ? TMath::Sqrt(errMx[5]) : -13,
		    dca.charge()*dca.pt(),    (errMx[9] >= 0 && dca.pt() > 0)  ? 100*TMath::Sqrt(errMx[9])*dca.pt() : -13,
		    dca.tanDip(),    (errMx[14] >= 0) ? TMath::Sqrt(errMx[14]): -13);
}

void   StDcaGeometry::Print(Option_t *option) const {cout << *this << endl;}

