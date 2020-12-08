# ==========================================
# CUSTUMIZED NUMERICAL RECICPES LIBRARY make
# ==========================================
# Compilation incompatibilities in functions:
# void newt(float x[], int n, int *check, void (*vecfunc)(int, float [], float []));
# void broydn(float x[], int n, int *check, void (*vecfunc)(int, float [], float []));
# Only true functions are compiled. Files containing main functions for examples are
# not compiled:
# rlftfrag.c  xran.c  fredex.c  xfit.c  badluk.c  sphoot.c  sphfpt.c  c_standard.c
# sfroid.c
# ==========================================
# Gnu C make to build Numerical Recipes Library (plus some personal additions)
# The generated libraries are libNum_Rec.a and libNRUTL.a
# ==========================================
# IMPORTANT NOTE:
# You need to compile and build up two Numerical Recipes libraries. Please do it 
# in two steps:
# . ~$ make 
# . ~$ make nrutil
#
PROJECT_DIR=${HOME}/PROJECT_STOCHASTIC_DIFFUSION/
C_NUMREC="$(PROJECT_DIR)/Library/C_NumRecip/"
EXTRA_INCLUDE_PATH= -I${C_NUMREC} 
CFLAGS=-g  -fPIE -DHAVE_INLINE ${EXTRA_INCLUDE_PATH}
###CFLAGS= -Wall -O21
COMPILER=gcc

OBJECTS = addint.o airy.o amebsa.o amoeba.o amotry.o amotsa.o anneal.o anorm2.o arcmak.o arcode.o arcsum.o asolve.o atimes.o avevar.o balanc.o banbks.o bandec.o banmul.o bcucof.o bcuint.o beschb.o bessi0.o bessi1.o bessi.o bessik.o bessj0.o bessj1.o bessj.o bessjy.o bessk0.o bessk1.o bessk.o bessy0.o bessy1.o bessy.o beta.o betacf.o betai.o bico.o bksub.o bnldev.o brent.o bsstep.o caldat.o chder.o chebev.o chebft.o chebpc.o chint.o chixy.o choldc.o cholsl.o chsone.o chstwo.o cisi.o cntab1.o cntab2.o convlv.o copy.o correl.o cosft1.o cosft2.o covsrt.o crank.o cyclic.o daub4.o dawson.o dbrent.o ddpoly.o decchk.o df1dim.o dfour1.o dfpmin.o dfridr.o dftcor.o dftint.o difeq.o dlinmin.o dpythag.o drealft.o dsprsax.o dsprstx.o dsvbksb.o dsvdcmp.o eclass.o eclazz.o ei.o eigsrt.o elle.o ellf.o ellpi.o elmhes.o erfcc.o erff.o erffc.o eulsum.o evlmem.o expdev.o expint.o f1dim.o factln.o factrl.o fasper.o fdjac.o fgauss.o fill0.o fit.o fitexy.o fixrts.o fleg.o flmoon.o fmin.o four1.o fourew.o fourfs.o fourn.o fpoly.o fred2.o fredin.o frenel.o frprmn.o ftest.o gamdev.o gamma_Double.o gammln.o gammp.o gammq.o gasdev.o gaucof.o gauher.o gaujac.o gaulag.o gauleg.o gaussj.o gcf.o gdf_double.o golden.o gser.o gser_double.o hpsel.o hpsort.o hqr.o hufapp.o hufdec.o hufenc.o hufmak.o hunt.o hypdrv.o hypgeo.o hypser.o icrc1.o icrc.o igray.o iindexx.o indexx.o interp.o irbit1.o irbit2.o jacobi.o jacobn.o julday.o kendl1.o kendl2.o kermom.o ks2d1s.o ks2d2s.o ksone.o kstwo.o laguer.o lfit.o linbcg.o linmin.o lnsrch.o locate.o lop.o lubksb.o ludcmp.o machar.o matadd.o matsub.o medfit.o memcof.o metrop.o mgfas.o mglin.o midexp.o midinf.o midpnt.o midsql.o midsqu.o miser.o mmid.o mnbrak.o mnewt.o moment.o mp2dfr.o mpdiv.o mpinv.o mpmul.o mpops.o mppi.o mprove.o mpsqrt.o mrqcof.o mrqmin.o odeint.o orthog.o pade.o pccheb.o pcshft.o pearsn.o period.o piksr2.o piksrt.o pinvs.o plgndr.o poidev.o polcoe.o polcof.o poldiv.o polin2.o polint.o powell.o predic.o probks.o psdes.o pwt.o pwtset.o pythag.o pzextr.o qgaus.o qrdcmp.o qromb.o qromo.o qroot.o qrsolv.o qrupdt.o qsimp.o qtrap.o quad3d.o quadct.o quadmx.o quadvl.o ran0.o ran1.o ran2.o ran3.o ran4.o rank.o ranpt.o ratint.o ratlsq.o ratval.o rc.o rd.o realft.o rebin.o red.o relax2.o relax.o resid.o revcst.o reverse.o rf.o rj.o rk4.o rkck.o rkdumb.o rkqs.o rlft3.o rofunc.o rotate.o rsolv.o rstrct.o rtbis.o rtflsp.o rtnewt.o rtsafe.o rtsec.o rzextr.o savgol.o scrsho.o select_nr.o selip.o shell.o shoot.o shootf.o simp1.o simp2.o simp3.o simplx.o simpr.o sinft.o slvsm2.o slvsml.o sncndn.o snrm.o sobseq.o solvde.o sor.o sort2.o sort3.o sort.o spctrm.o spear.o spectral_Density.o sphbes.o splie2.o splin2.o spline.o splint.o spread.o sprsax.o sprsin.o sprspm.o sprstm.o sprstp.o sprstx.o stifbs.o stiff.o stoerm.o svbksb.o svdcmp.o svdfit.o svdvar.o toeplz.o tptest.o tqli.o trapzd.o tred2.o tridag.o tridag_Function.o trncst.o trnspt.o ttest.o tutest.o twofft.o vander.o vegas.o voltra.o wt1.o wtn.o wwghts.o zbrac.o zbrak.o zbrent.o zrhqr.o zriddr.o zroots.o

.c.o:
	${COMPILER} $(CFLAGS) -c $<

NumRec: ${OBJECTS} nr.h
	ar -q libNumRec.a ${OBJECTS}

nrutil: nrutil.o nrutil.h
	ar -q libNRUTL.a nrutil.o

all: NumRec nrutil

clean:  
	rm -f *.o *~ 
