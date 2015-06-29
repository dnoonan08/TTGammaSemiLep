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
	if 'nJets' in varname:
		sfname = 'ttgamma fraction'
	signalFractionVar = RooRealVar(sfname,sfname, 0.8, 0.0, 1.0)
	sumPdf = RooAddPdf('totalPdf','signal and background', signalPdf, backgroundPdf, signalFractionVar)
	
	# fit
	sumPdf.fitTo( dataDataHist, RooFit.SumW2Error(kFALSE), RooFit.PrintLevel(-1) )
	
	if plotName!='':
		# plot results
		c1 = TCanvas('c1', 'c1', 800, 600)
		plotter = RooPlot('myplot','',sihihVar,varmin,varmax,20) # nBins is dummy
		dataDataHist.plotOn(plotter, RooFit.Name('data'))
		sumPdf.plotOn(plotter, RooFit.Name('sum'), RooFit.LineColor(kRed))
		sumPdf.plotOn(plotter,RooFit.Components('signalPdf'), RooFit.Name('signal'), 
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
        qcdVar = RooRealVar(qcdfname,qcdfname, qcdIntegral,0.5*qcdIntegral,1.5*qcdIntegral)
        otherMCVar = RooRealVar(otherMCfname, otherMCfname,otherMCIntegral,0.8*otherMCIntegral,1.2*otherMCIntegral) 


	Gauss_QCD = RooGaussian("gauss_QCD","gauss_QCD",qcdVar,RooFit.RooConst(qcdIntegral),RooFit.RooConst(0.5*qcdIntegral))
        Gauss_otherMC =  RooGaussian("gauss_otherMC","gauss_otherMC",otherMCVar,RooFit.RooConst(otherMCIntegral),RooFit.RooConst(.2*otherMCIntegral))

        #qcdVar.setConstant(True)
        #otherMCVar.setConstant(True)
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
        sumPdf.fitTo( dataDataHist, RooFit.Extended(True), RooFit.SumW2Error(kFALSE), RooFit.PrintLevel(-1) )

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
        print 'fit returned value ',signalVar.getVal(),' +- ',signalVar.getError()
        return (signalVar.getVal(),signalVar.getError(),   bkgVar.getVal(),bkgVar.getError(), otherMCVar.getVal(),otherMCVar.getError(), qcdVar.getVal(), qcdVar.getError())

                                                                            
openfiles = {}

def get1DHist(filename, histname):
	if filename not in openfiles:
		openfiles[filename] =TFile(filename,'READ')
	file = openfiles[filename]
		
	hist = file.Get(histname)
	hist.SetDirectory(0)
	hist.SetFillColor(0)
	#hist.Sumw2()
	return hist


M3file_photon = 'templates_barrel_scaled.root'
M3file_presel_scaled = 'templates_presel_scaled.root'


def integral_bins(hist,bin1,bin2):
	err = Double(0.0)
	integr = hist.IntegralAndError(bin1, bin2, err)
	print integr,err
	return integr,err
def doNjetsfit_photon():
	print 'nJets fit after photon selection'
	varToFit = 'nJets'
	lowfitBin = 4
	highfitBin = 9
	
	DataHist = get1DHist(M3file_photon, 'Data_'+varToFit)
	
	TopHist = get1DHist(M3file_photon, 'TTGamma_'+varToFit)
	
	TTJetsHist = get1DHist(M3file_photon, 'TTJets_'+varToFit)
	#TopHist.Add(TTJetsHist)
	
	ZJetsHist = get1DHist(M3file_photon, 'ZJets_'+varToFit)
	#ZJetsHist = get1DHist(M3file_photon, 'ZJets_signal_'+varToFit)
	#ZJetsHist.Add(get1DHist(M3file_photon, 'ZJets_fake_'+varToFit))
	#ZJetsHist.Add(get1DHist(M3file_photon, 'ZJets_electron_'+varToFit).Scale(1.5))
	
	WJetsHist = get1DHist(M3file_photon, 'WJets_'+varToFit)
	SingleTopHist = get1DHist(M3file_photon, 'SingleTop_'+varToFit)
	QCDHist = get1DHist(M3file_photon, 'QCD_'+varToFit)
	
	# remove extra contributions from data
	DataHist.Add(TTJetsHist, -1.0)
	DataHist.Add(WJetsHist, -1.0)
	DataHist.Add(ZJetsHist, -1.0)
	DataHist.Add(SingleTopHist, -1.0)
	DataHist.Add(QCDHist, -1.0)
	
	BGHist = get1DHist(M3file_photon, 'Vgamma_'+varToFit)
	
	dataInt,dataIntErr = integral_bins(DataHist,lowfitBin,highfitBin)
	topInt,topIntErr = integral_bins(TopHist,lowfitBin,highfitBin)
	bgInt,bgIntErr = integral_bins(BGHist,lowfitBin,highfitBin)
	
	(nJetsTopFrac, nJetsTopFracErr) = makeFit(varToFit+', photon selection', 2.5, 8.5, TopHist, BGHist, DataHist, 'plots/'+varToFit+'_photon_fit.png')

	bgSF = (1.0 - nJetsTopFrac) * dataInt / bgInt
	bgSFerror = bgSF * ( (nJetsTopFracErr / (1.0 - nJetsTopFrac))**2 + (bgIntErr/bgInt)**2 + (dataIntErr/dataInt)**2 )**0.5
	
	topSF = nJetsTopFrac * dataInt / topInt
	topSFErr = topSF * ( (nJetsTopFracErr / nJetsTopFrac)**2 + (topIntErr/topInt)**2 + (dataIntErr/dataInt)**2 )**0.5
	
	print '#'*80
	print 'Correction to Vgamma samples after nJets fit: ', bgSF, ' +-', bgSFerror, '(fit + stat error)'
	print 'Correction to ttgamma sample after nJets fit: ', topSF, ' +-', topSFErr, '(fit + stat error)'
	print '#'*80
	return (bgSF,bgSFerror)



def doM3fit_photon():
	print 'M3 fit after photon selection'
	varToFit = 'M3'
	DataHist = get1DHist(M3file_photon, 'Data_'+varToFit)

	
	# ttjets and ttgamma shapes after photon selection
	TopHist = get1DHist(M3file_photon, 'TTJets_'+varToFit)
	#TopHist = get1DHist(M3file_presel_scaled, 'TTJets_'+varToFit)
	TopHist.Add(get1DHist(M3file_photon, 'TTGamma_'+varToFit))
	WJetsHist = get1DHist(M3file_presel_scaled, 'WJets_'+varToFit)
	WJetsHist.Scale(get1DHist(M3file_photon, 'WJets_'+varToFit).Integral() / WJetsHist.Integral())
		
	ZJetsHist = get1DHist(M3file_presel_scaled, 'ZJets_'+varToFit)
	ZJetsHist.Scale(get1DHist(M3file_photon, 'ZJets_'+varToFit).Integral() / ZJetsHist.Integral())
		
	SingleTopHist = get1DHist(M3file_presel_scaled, 'SingleTop_'+varToFit)
	SingleTopHist.Scale(get1DHist(M3file_photon, 'SingleTop_'+varToFit).Integral() / SingleTopHist.Integral())
	print M3file_presel_scaled
	QCDHist = get1DHist(M3file_presel_scaled, 'QCD_'+varToFit)
	QCDHist.Scale(get1DHist(M3file_photon, 'QCD_'+varToFit).Integral() / QCDHist.Integral())
	
	
	VgHist = get1DHist(M3file_presel_scaled, 'Vgamma_'+varToFit)
	VgHist.Scale(get1DHist(M3file_photon, 'Vgamma_'+varToFit).Integral() / VgHist.Integral())
	
	# add all backgrounds
	otherMCHist = ZJetsHist
	otherMCHist.Add(SingleTopHist)
	otherMCHist.Add(VgHist)

	
	DataHist.Rebin(2)
	TopHist.Rebin(2)
	WJetsHist.Rebin(2)
        otherMCHist.Rebin(2)
	QCDHist.Rebin(2)

	(m3Top, m3TopErr,m3Wjets, m3WjetsErr, m3otherMC, m3otherMCerr, m3QCD, m3QCDerr) = makenewFit(varToFit+'(GeV), photon selection', 0.0,800.0, TopHist, WJetsHist,otherMCHist, QCDHist, DataHist, 'plots/'+varToFit+'_photon_fit.png')
        lowfitBin = DataHist.FindBin(0.01)
        highfitBin = DataHist.FindBin(799.99)

        dataInt = DataHist.Integral(lowfitBin,highfitBin)
        topInt = TopHist.Integral(lowfitBin,highfitBin)
        WJInt = WJetsHist.Integral(lowfitBin,highfitBin)
        QCDInt = QCDHist.Integral(lowfitBin,highfitBin)
	otherMCInt = otherMCHist.Integral(lowfitBin, highfitBin)
        TopSF = m3Top/ topInt
        TopSFerror = m3TopErr/ topInt


	print
	print '#'*80
	print 'Total amount of Top events in fit:', m3Top, '+-', m3TopErr
	print 'Total amount of WJets events in fit:', m3Wjets, '+-', m3WjetsErr
	print 'Total amount of QCD events in fit:', m3QCD, '+-', m3QCDerr
	print 'Total amount of Other MC events in fit:', m3otherMC, '+-', m3otherMCerr
	print '#'*80

	print


        print '#'*80
        print 'Correction to the Top scale factor: ', TopSF, ' +-', TopSFerror, '(fit error only)'
        WJetsSF = m3Wjets / WJInt
        WJetsSFerror = m3WjetsErr / WJInt
	otherMCSF = m3otherMC/otherMCInt
	otherMCSFerror = m3otherMCerr/otherMCInt
        #QCDSF = (1-m3Wjets-m3TopFrac)*dataInt/QCDInt
	totMC = m3Top + m3Wjets + m3otherMC + m3QCD
	m3_topFrac = m3Top/totMC
	m3_topFracErr = m3TopErr/totMC	
	print 'the top fraction after photon selection :',m3_topFrac
	print 'the top fraction error after photon selection:', m3_topFracErr 	
	print 'Correction to WJets scale factor: ', WJetsSF, ' +-',WJetsSFerror,'(fit error only)'
	print 'Correction to otherMC scale factor: ', otherMCSF, ' +-',otherMCSFerror,'(fit error only)'
        print '#'*80
	######## Calculate the top fraction (topEvents over total MC events) and return it as well ########
        return (TopSF, TopSFerror, WJetsSF, WJetsSFerror,otherMCSF, otherMCSFerror, m3_topFrac, m3_topFracErr) 

	
def doQCDfit_photon():
        varToFit = 'MET'
        DataHist = get1DHist(M3file_photon, 'Data_'+varToFit)


        # ttjets and ttgamma shapes after photon selection
        TopHist = get1DHist(M3file_photon, 'TTJets_'+varToFit)
        #TopHist = get1DHist(M3file_presel_scaled, 'TTJets_'+varToFit)
        TopHist.Add(get1DHist(M3file_photon, 'TTGamma_'+varToFit))


        # other shapes from top selection but scaled according to event count in photon selection
        WJetsHist = get1DHist(M3file_presel_scaled, 'WJets_'+varToFit)
        WJetsHist.Scale(get1DHist(M3file_photon, 'WJets_'+varToFit).Integral() / WJetsHist.Integral())

        ZJetsHist = get1DHist(M3file_presel_scaled, 'ZJets_'+varToFit)
        ZJetsHist.Scale(get1DHist(M3file_photon, 'ZJets_'+varToFit).Integral() / ZJetsHist.Integral())

        SingleTopHist = get1DHist(M3file_presel_scaled, 'SingleTop_'+varToFit)
        SingleTopHist.Scale(get1DHist(M3file_photon, 'SingleTop_'+varToFit).Integral() / SingleTopHist.Integral())
        QCDHist = get1DHist(M3file_presel_scaled, 'QCD_'+varToFit)
        QCDHist.Scale(get1DHist(M3file_photon, 'QCD_'+varToFit).Integral() / QCDHist.Integral())


        VgHist = get1DHist(M3file_presel_scaled, 'Vgamma_'+varToFit)
        VgHist.Scale(get1DHist(M3file_photon, 'Vgamma_'+varToFit).Integral() / VgHist.Integral())
	# add all backgrounds
        MCHist = WJetsHist
        #BGHist = SingleTopHist
        MCHist.Add(ZJetsHist)
        MCHist.Add(TopHist)
        MCHist.Add(SingleTopHist)
        MCHist.Add(VgHist)



        (metFrac, metFracErr) = makeFit(varToFit, 0., 200.0, QCDHist, MCHist, DataHist, 'plots/'+varToFit+'_QCD_fit_photon.png')
        # recalculate data-driven QCD normalization
        lowbin = DataHist.FindBin(0.01)
        highbin = DataHist.GetNbinsX()+1# overflow bin included
      
        print 'Will calculate integral in the bin range:', lowbin, highbin
        dataInt = DataHist.Integral(lowbin, highbin)
        dataIntwithMETcut = DataHist.Integral((DataHist.FindBin(20.01)),highbin)
        qcdIntwithMETcut = QCDHist.Integral((DataHist.FindBin(20.01)),highbin)
        print 'Integral of Data total: ', dataInt
        qcdInt = QCDHist.Integral(lowbin, highbin)
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
        QCDSF_photon = metFrac*dataInt/qcdInt
        QCDSFerror_photon = metFracErr*dataInt/qcdInt
        print 'Amount of QCD in signal :' ,qcdIntwithMETcut*QCDSF_photon,'+-',qcdIntwithMETcut*QCDSFerror_photon

        print 'Scale factor for QCD in nominal MET range: ', QCDSF_photon,' +-',QCDSFerror_photon,'(fit error only)'
        print 'Correction to all MC scale factors: ', (1-metFrac)*dataInt/mcInt, ' +-',metFracErr*dataInt/mcInt,'(fit error only)'
        print '#'*80
        return (QCDSF_photon, QCDSFerror_photon)

def doQCDlowfit_photon():
        varToFit = 'MET_low'
        DataHist = get1DHist(M3file_photon, 'Data_'+varToFit)


        # ttjets and ttgamma shapes after photon selection
        TopHist = get1DHist(M3file_photon, 'TTJets_'+varToFit)
        #TopHist = get1DHist(M3file_presel_scaled, 'TTJets_'+varToFit)
        TopHist.Add(get1DHist(M3file_photon, 'TTGamma_'+varToFit))


        # other shapes from top selection but scaled according to event count in photon selection
        WJetsHist = get1DHist(M3file_presel_scaled, 'WJets_'+varToFit)
        WJetsHist.Scale(get1DHist(M3file_photon, 'WJets_'+varToFit).Integral() / WJetsHist.Integral())

        ZJetsHist = get1DHist(M3file_presel_scaled, 'ZJets_'+varToFit)
        ZJetsHist.Scale(get1DHist(M3file_photon, 'ZJets_'+varToFit).Integral() / ZJetsHist.Integral())

        SingleTopHist = get1DHist(M3file_presel_scaled, 'SingleTop_'+varToFit)
        SingleTopHist.Scale(get1DHist(M3file_photon, 'SingleTop_'+varToFit).Integral() / SingleTopHist.Integral())
        QCDHist = get1DHist(M3file_presel_scaled, 'QCD_'+varToFit)
        QCDHist.Scale(get1DHist(M3file_photon, 'QCD_'+varToFit).Integral() / QCDHist.Integral())


        VgHist = get1DHist(M3file_presel_scaled, 'Vgamma_'+varToFit)
        VgHist.Scale(get1DHist(M3file_photon, 'Vgamma_'+varToFit).Integral() / VgHist.Integral())
        # add all backgrounds
        MCHist = WJetsHist
        #BGHist = SingleTopHist
        MCHist.Add(ZJetsHist)
        MCHist.Add(TopHist)
        MCHist.Add(SingleTopHist)
        MCHist.Add(VgHist)
        (metFrac, metFracErr) = makeFit(varToFit, 0., 20.0, QCDHist, MCHist, DataHist, 'plots/'+varToFit+'_QCD_fit_photon.png')
        # recalculate data-driven QCD normalization
        lowbin = DataHist.FindBin(0.01)
        highbin = DataHist.FindBin(19.99)# overflow bin included
        print 'Will calculate integral in the bin range:', lowbin, highbin
	print 'MET fit in 0->20 GeV range'
        dataInt = DataHist.Integral(lowbin, highbin)
        print 'Integral of Data total: ', dataInt
        qcdInt = QCDHist.Integral(lowbin, highbin)
        mcInt = MCHist.Integral(lowbin, highbin)
        print 'Integral of data-driven QCD in the desired range: ', qcdInt
        print '#'*80
        # take into account only fit error
        # stat errors on histograms are treated while calculating the final answer
        print 'the met fraction:', metFrac
        QCDamount = metFrac*dataInt
	QCDerror = metFracErr*dataInt
        print 'Total Amount of QCD :',QCDamount,'+-',QCDerror
        QCD_low_SF_photon = metFrac*dataInt/qcdInt
        QCD_low_SFerror_photon = metFracErr*dataInt/qcdInt

        print 'Scale factor for QCD in nominal MET range: ', QCD_low_SF_photon,' +-',QCD_low_SFerror_photon,'(fit error only)'
        print 'Correction to all MC scale factors: ', (1-metFrac)*dataInt/mcInt, ' +-',metFracErr*dataInt/mcInt,'(fit error only)'
        print '#'*80
        return (QCD_low_SF_photon, QCD_low_SFerror_photon)
                                 
