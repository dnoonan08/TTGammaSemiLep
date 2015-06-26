#! /usr/bin/env python
from ROOT import *
from ROOT import RooFit
import sys
import os
import math
import re

def makeFit(varname, varmin, varmax, signalHist, backgroundHist, dataHist, plotName):
	# RooFit variables
	sihihVar = RooRealVar(varname, varname, varmin, varmax)
	sihihArgList = RooArgList()
	sihihArgList.add(sihihVar)
	sihihArgSet = RooArgSet()
	sihihArgSet.add(sihihVar)

	# create PDFs
	signalDataHist = RooDataHist('signalDataHist','signal RooDataHist', sihihArgList, signalHist)
	signalPdf = RooHistPdf('signalPdf',varname+' of signal', sihihArgSet, signalDataHist)

	backgroundDataHist = RooDataHist('backgroundDataHist','background RooDataHist', sihihArgList, backgroundHist)
	backgroundPdf = RooHistPdf('backgroundPdf',varname+' of background', sihihArgSet, backgroundDataHist)



	# data
	dataDataHist = RooDataHist('data '+varname, varname+' in Data', sihihArgList, dataHist)

	# signal fraction parameter
	sfname = 'signal fraction'
	if 'MET' in varname:
		sfname = 'multijet fraction'
	if 'M3' in varname:
		sfname = 'ttbar fraction'
	signalFractionVar = RooRealVar(sfname,sfname, 0.5, 0.0, 1.0)
	sumPdf = RooAddPdf('totalPdf','signal+background', signalPdf, backgroundPdf, signalFractionVar)
	
	# fit
	sumPdf.fitTo( dataDataHist, RooFit.SumW2Error(0), RooFit.PrintLevel(-1) )
	
	if plotName!='':
		# plot results
		c1 = TCanvas('c1', 'c1', 800, 600)
		plotter = RooPlot('myplot','',sihihVar,varmin,varmax,20) # nBins is dummy
		dataDataHist.plotOn(plotter, RooFit.Name('data'))
		sumPdf.plotOn(plotter, RooFit.Name('sum'), RooFit.LineColor(kRed))
		sumPdf.plotOn(plotter, RooFit.Components('signalPdf'), RooFit.Name('signal'), 
			RooFit.LineColor(50))
		sumPdf.plotOn(plotter, RooFit.Components('backgroundPdf'), RooFit.Name('background'), 
			RooFit.LineColor(kBlue))
		sumPdf.paramOn(plotter) # fix

		plotter.Draw()
		plotter.GetYaxis().SetTitleOffset(1.4)
		plotter.GetXaxis().SetTitle("GeV")
		c1.SaveAs(plotName)
        
	print 'fit returned value ',signalFractionVar.getVal(),' +- ',signalFractionVar.getError()
	return (signalFractionVar.getVal(),signalFractionVar.getError())


def makenewFit(varname, varmin, varmax, signalHist, backgroundHist, otherMCHist, qcdHist, dataHist, plotName):
        # RooFit variables
        sihihVar = RooRealVar(varname, varname, varmin, varmax)
        sihihArgList = RooArgList()
        sihihArgList.add(sihihVar)
        sihihArgSet = RooArgSet()
        sihihArgSet.add(sihihVar)

        # create PDFs
        signalDataHist = RooDataHist('signalDataHist','signal RooDataHist', sihihArgList, signalHist)
        signalPdf = RooHistPdf('signalPdf',varname+' of signal', sihihArgSet, signalDataHist)

        backgroundDataHist = RooDataHist('backgroundDataHist','background RooDataHist', sihihArgList, backgroundHist)
        bkgPdf = RooHistPdf('backgroundPdf',varname+' of background', sihihArgSet, backgroundDataHist)
       
        qcdDataHist = RooDataHist('qcdDataHist','qcd RooDataHist', sihihArgList, qcdHist)
        qcdPdf = RooHistPdf('qcdPdf',varname+' of qcd', sihihArgSet, qcdDataHist)
    
        otherMCDataHist = RooDataHist('otherMCDataHist','otherMC RooDataHist', sihihArgList, otherMCHist)
        otherMCPdf = RooHistPdf('otherMCPdf',varname+' of otherMC', sihihArgSet, otherMCDataHist)


 
        

        # data
        dataDataHist = RooDataHist('data '+varname, varname+' in Data', sihihArgList, dataHist)

        # signal fraction parameter
        sfname = 'signal fraction'
        if 'MET' in varname:
                sfname = 'multijet fraction'
        if 'M3' in varname:
                sfname = 'ttbar total'
                bkgfname = 'background total'
                qcdfname = ' qcd total'
		otherMCfname = 'other MC total'
	signalIntegral   = signalHist.Integral()
	bkgIntegral      = backgroundHist.Integral()
	qcdIntegral      = qcdHist.Integral()
	otherMCIntegral  = otherMCHist.Integral()
        signalVar = RooRealVar(sfname,sfname, signalIntegral,0.,5.*signalIntegral)
        bkgVar = RooRealVar(bkgfname,bkgfname, bkgIntegral,0.,5.*bkgIntegral)
        qcdVar = RooRealVar(qcdfname,qcdfname, qcdIntegral,0.*qcdIntegral,5.*qcdIntegral)
        otherMCVar = RooRealVar(otherMCfname, otherMCfname,otherMCIntegral,0.*otherMCIntegral,5.*otherMCIntegral) 
	#constraints:
        #qcdVar.setVal(qcd)
        qcdVar.setConstant(True)

	#otherMCVar.setVal(otherMC)
	#otherMCVar.setConstant(True)
        #setting QCD and otherMC with a gaussian constraint:
        Gauss_QCD = RooGaussian("gauss_QCD","gauss_QCD",qcdVar,RooFit.RooConst(qcdIntegral),RooFit.RooConst(0.5*qcdIntegral))
	Gauss_otherMC =  RooGaussian("gauss_otherMC","gauss_otherMC",otherMCVar,RooFit.RooConst(otherMCIntegral),RooFit.RooConst(.2*otherMCIntegral))

	constraints = RooArgSet(Gauss_QCD,Gauss_otherMC)
        listPdfs = RooArgList(signalPdf,\
                                   bkgPdf,\
                                   otherMCPdf,\
                                   qcdPdf)
        listcoeff = RooArgList(signalVar,\
                                  bkgVar,\
                                  otherMCVar,\
				  qcdVar)


	sumPdf = RooAddPdf('totalPdf','signal+background+qcd', listPdfs, listcoeff )
        
        # fit
        sumPdf.fitTo( dataDataHist, RooFit.InitialHesse(True),RooFit.Minos(True),RooFit.ExternalConstraints(constraints), RooFit.Extended(True), RooFit.SumW2Error(kFALSE), RooFit.PrintLevel(-1) )

        if plotName!='':
                # plot results
                c1 = TCanvas('c1', 'c1', 800, 600)
                plotter = RooPlot('myplot','',sihihVar,varmin,varmax,20) # nBins is dummy
                dataDataHist.plotOn(plotter, RooFit.Name('data'))
                sumPdf.plotOn(plotter, RooFit.Name('sum'), RooFit.LineColor(kRed))
                sumPdf.plotOn(plotter, RooFit.Components('signalPdf'), RooFit.Name('signal'),
                    RooFit.LineColor(kGreen))
                sumPdf.plotOn(plotter, RooFit.Components('backgroundPdf'), RooFit.Name('background'),
                    RooFit.LineColor(kBlue))
		sumPdf.plotOn(plotter, RooFit.Components('otherMCPdf'), RooFit.Name('otherMC'),
                    RooFit.LineColor(6))

                sumPdf.plotOn(plotter, RooFit.Components('qcdPdf'), RooFit.Name('qcd'),
                    RooFit.LineColor(50))

                sumPdf.paramOn(plotter) # fix

                plotter.Draw()
                plotter.GetYaxis().SetTitleOffset(1.4)
		plotter.GetXaxis().SetTitle("GeV")
                c1.SaveAs(plotName)
        print 'fit returned value signal:',signalVar.getVal(),' +- ',signalVar.getError()
        
        return (signalVar.getVal(),signalVar.getError(),   bkgVar.getVal(),bkgVar.getError(), otherMCVar.getVal(),otherMCVar.getError(),qcdVar.getVal(),qcdVar.getError())


                                                                                                                                                      
                    
openfiles = {}

def get1DHist(filename, histname):
	if filename not in openfiles:
		openfiles[filename] = TFile(filename,'READ')
	file = openfiles[filename]
		
	hist = file.Get(histname)
	hist.SetDirectory(0)
	hist.SetFillColor(0)
	#hist.Sumw2()
	return hist

qcdMETfile = 'templates_presel_nomet_qcd.root'
normMETfile = 'templates_presel_nomet.root'

M3file = 'templates_presel.root'
#M3file = 'templates_presel_nomet.root'
M3file_photon = 'templates_barrel_scaled.root'
M3file_presel_scaled = 'templates_presel_scaled.root'
def doQCD_lowfit():
	varToFit = 'MET_low'
	qcdDataHist = get1DHist(qcdMETfile, 'Data_'+varToFit)
        # remove MC contribution
        qcdDataHist.Add(get1DHist(qcdMETfile, 'TTJets_'+varToFit), -1)
        qcdDataHist.Add(get1DHist(qcdMETfile, 'WJets_'+varToFit), -1)
        qcdDataHist.Add(get1DHist(qcdMETfile, 'SingleTop_'+varToFit), -1)

        DataHist = get1DHist(normMETfile, 'Data_'+varToFit)

        MCHist = get1DHist(normMETfile, 'TTJets_'+varToFit)
        MCHist.Add(get1DHist(normMETfile, 'TTGamma_'+varToFit))
        MCHist.Add(get1DHist(normMETfile, 'WJets_'+varToFit))
        MCHist.Add(get1DHist(normMETfile, 'ZJets_'+varToFit))
        MCHist.Add(get1DHist(normMETfile, 'SingleTop_'+varToFit))

        (metFrac, metFracErr) = makeFit(varToFit, 0., 20.0, qcdDataHist, MCHist, DataHist, 'plots/'+varToFit+'_QCD_fit.png')
        # recalculate data-driven QCD normalization
        lowbin = DataHist.FindBin(0.01)
        highbin =  DataHist.FindBin(19.99)# overflow bin included
        print 'Will calculate integral in the bin range:', lowbin, highbin
	print 'MET fit in region 0 ->20 GeV after photon selection:'

        dataInt = DataHist.Integral(lowbin, highbin)
        print 'Integral of Data total in low MET region: ', dataInt
        qcdInt = qcdDataHist.Integral(lowbin, highbin)
        mcInt = MCHist.Integral(lowbin, highbin)
        print 'Integral of data-driven QCD low MET range: ', qcdInt
        print '#'*80
        print 'the met fraction:', metFrac
        QCDamount = metFrac*dataInt
        QCDerror = metFracErr*dataInt
        print 'Total Amount of QCD :',QCDamount,'+-',QCDerror
        QCD_low_SF = metFrac*dataInt/qcdInt
        QCD_low_SFerror = metFracErr*dataInt/qcdInt

     
        print 'Scale factor for QCD in nominal MET range: ', QCD_low_SF,' +-',QCD_low_SFerror,'(fit error only)'
        print 'Correction to all MC scale factors: ', (1-metFrac)*dataInt/mcInt, ' +-',metFracErr*dataInt/mcInt,'(fit error only)'
        print '#'*80
	return (QCD_low_SF, QCD_low_SFerror)

def doQCDfit():
	varToFit = 'MET'

	qcdDataHist = get1DHist(qcdMETfile, 'Data_'+varToFit)
	# remove MC contribution
	qcdDataHist.Add(get1DHist(qcdMETfile, 'TTJets_'+varToFit), -1)
	qcdDataHist.Add(get1DHist(qcdMETfile, 'WJets_'+varToFit), -1)
        qcdDataHist.Add(get1DHist(qcdMETfile, 'SingleTop_'+varToFit), -1)

	DataHist = get1DHist(normMETfile, 'Data_'+varToFit)

	MCHist = get1DHist(normMETfile, 'TTJets_'+varToFit)
	MCHist.Add(get1DHist(normMETfile, 'TTGamma_'+varToFit))
	MCHist.Add(get1DHist(normMETfile, 'WJets_'+varToFit))
	MCHist.Add(get1DHist(normMETfile, 'ZJets_'+varToFit))
	MCHist.Add(get1DHist(normMETfile, 'SingleTop_'+varToFit))

	(metFrac, metFracErr) = makeFit(varToFit, 0., 200.0, qcdDataHist, MCHist, DataHist, 'plots/'+varToFit+'_QCD_fit.png')
	# recalculate data-driven QCD normalization
	lowbin = DataHist.FindBin(0.01)
	highbin =   DataHist.GetNbinsX()+1 # overflow bin included
        print DataHist.GetBinContent(3), DataHist.GetBinContent(4), DataHist.GetBinContent(5)
	print 'MET fit in region 0 ->200 GeV: after photon selection'
	print 'Will calculate integral in the bin range:', lowbin, highbin
	dataInt = DataHist.Integral(lowbin, highbin)
	dataIntwithMETcut = DataHist.Integral(DataHist.FindBin(20.01),highbin)
        qcdIntwithMETcut = qcdDataHist.Integral(DataHist.FindBin(20.01),highbin)
	print 'Integral of Data total: ', dataInt
	qcdInt = qcdDataHist.Integral(lowbin, highbin)
	mcInt = MCHist.Integral(lowbin, highbin)
	print 'Integral of data-driven QCD in the desired range: ', qcdInt
	print '#'*80
	# take into account only fit error
	# stat errors on histograms are treated while calculating the final answer
	print 'the met fraction:', metFrac
	print 'the amount of Data with MET cut:' ,dataIntwithMETcut
	QCDamount = metFrac*dataInt
	QCDerror = metFracErr*dataInt
	print 'Total Amount of QCD without MET cut:',QCDamount,'+-',QCDerror
	QCDSF = metFrac*dataInt/qcdInt
	QCDSFerror = metFracErr*dataInt/qcdInt

        print 'Amount of QCD in signal :' ,qcdIntwithMETcut*QCDSF,'+-',qcdIntwithMETcut*QCDSFerror

	print 'Scale factor for QCD in nominal MET range: ', QCDSF,' +-',QCDSFerror,'(fit error only)'
	print 'Correction to all MC scale factors: ', (1-metFrac)*dataInt/mcInt, ' +-',metFracErr*dataInt/mcInt,'(fit error only)'
	print '#'*80
	return (QCDSF, QCDSFerror)



def doM3fit():
	print 'now do M3 fit'

	varToFit = 'M3'

	DataHist = get1DHist(M3file, 'Data_'+varToFit)

	TopHist = get1DHist(M3file, 'TTJets_'+varToFit)
	TopHist.Add(get1DHist(M3file, 'TTGamma_'+varToFit))



	WJHist = get1DHist(M3file, 'WJets_'+varToFit)
	
        
	otherMCHist= get1DHist(M3file, 'ZJets_'+varToFit)
	otherMCHist.Add(get1DHist(M3file, 'SingleTop_'+varToFit))
	otherMCHist.Add(get1DHist(M3file, 'Vgamma_'+varToFit))
	
	QCDHist = get1DHist(M3file, 'QCD_'+varToFit)

	(m3Top, m3TopErr, m3Wjets, m3WjetsErr, m3otherMC, m3otherMCerr,m3QCD,m3QCDerr) = makenewFit(varToFit+'(GeV)', 0.0, 800.0, TopHist, WJHist, otherMCHist, QCDHist, DataHist, 'plots/'+varToFit+'_fit.png')
	lowfitBin = DataHist.FindBin(0.01)
	highfitBin = DataHist.FindBin(799.99)
			
	dataInt = DataHist.Integral(lowfitBin,highfitBin)
	topInt = TopHist.Integral(lowfitBin,highfitBin)
	WJInt = WJHist.Integral(lowfitBin,highfitBin)
        QCDInt = QCDHist.Integral(lowfitBin,highfitBin)	
	otherMCInt = otherMCHist.Integral(lowfitBin,highfitBin)
	TopSF = m3Top/ topInt
	TopSFerror = m3TopErr/ topInt
	
	print '#'*80
	print 'Correction to the Top scale factor: ', TopSF, ' +-', TopSFerror, '(fit error only)'
	WJetsSF = m3Wjets / WJInt
	WJetsSFerror = m3WjetsErr / WJInt
	otherMCSF = m3otherMC/otherMCInt 
	otherMCSFerror = m3otherMCerr/otherMCInt
	QCDSF = m3QCD/QCDInt
        QCDSFerror = m3QCDerr/QCDInt
	#QCDSF = (1-m3Wjets-m3TopFrac)*dataInt/QCDInt
        #QCDSFErr = (1-m3WjetsErr-m3TopFracErr)*dataInt/QCDInt 
	print 'Correction to WJets scale factor: ', WJetsSF, ' +-',WJetsSFerror,'(fit error only)'
	print 'Correction to QCD scale factor: ', QCDSF, ' +-',QCDSFerror,'(fit error only)'
	print 'Correction to OtherMC scale factor: ', otherMCSF, ' +-',otherMCSFerror,'(fit error only)'
	print '#'*80
	return (TopSF, TopSFerror, WJetsSF, WJetsSFerror, otherMCSF, otherMCSFerror, QCDSF, QCDSFerror)

