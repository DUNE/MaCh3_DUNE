#!/bin/bash

wRC="true"

#WRC
if [ "$wRC" = "true" ]; then
    outfolder="contours_datasensitivitycomparisons_wRC"
    datacontourfile="jointFit2015_wFGD2_dataFit_run17c_wRC_valorMargin_contours_disapp.root"
    dataappcontourfile="jointFit2015_wFGD2_dataFit_run17c_wRC_valorMargin_contours_app.root"
    sensitivitycontourfile="contours_wRC_fixed1D_forofficial/contours_th23dm23_both.root"
    sensitivityappcontourfile="contours_wRC_fixed1D_forofficial/contours_th13dcp_both.root"
    datainputfile="jointFit2015_wFGD2_dataFit_run17c_wRC.root"
    sensitivityinputfile="jointFit2015_wFGD2_asimov1_run17c_wRC.root"
    leg1="Data (T2K only)"
    leg2="Sensitivity (T2K only)"
else
#WORC
    outfolder="contours_datasensitivitycomparisons_woRC"
    datacontourfile="jointFit2015_wFGD2_dataFit_run17c_woRC_SkdetBugfixed_valorMargin_contours_disapp.root"
    dataappcontourfile="jointFit2015_wFGD2_dataFit_run17c_woRC_SkdetBugfixed_valorMargin_contours_app.root"
    sensitivitycontourfile="contours_woRC_fixed1D_forofficial/contours_th23dm23_both.root"
    sensitivityappcontourfile="contours_woRC_fixed1D_forofficial/contours_th13dcp_both.root"
    datainputfile="jointFit2015_wFGD2_dataFit_run17c_woRC_SkdetBugfixed.root"
    sensitivityinputfile="jointFit2015_wFGD2_asimov1_run17c_woRC.root"
    legadd="(T2K+Reactor)"
fi


mkdir $outfolder

root -l -q 'utils/CompareContours.C++("'$datacontourfile'","'$sensitivitycontourfile'","Data '$legadd'","Sensitivity '$legadd'",false,"comparedcontours_th23dm23_datavssensitivity")'

mv comparedcontours_th23dm23_datavssensitivity.pdf $outfolder/comparedcontours_th23dm23_datavssensitivity.pdf
mv comparedcontours_th23dm23_datavssensitivity.png $outfolder/comparedcontours_th23dm23_datavssensitivity.png
mv comparedcontours_th23dm23_datavssensitivity.root $outfolder/comparedcontours_th23dm23_datavssensitivity.root

root -l -q -b 'utils/macros/makeOfficialPlotStyle.C++("'$outfolder'/comparedcontours_th23dm23_datavssensitivity.root","'$outfolder'/comparedcontours_th23dm23_datavssensitivity_official.root","T2K Run 1-7c preliminary","c1")'

root -l -q 'utils/CompareContours.C++("'$dataappcontourfile'","'$sensitivityappcontourfile'","Data '$legadd'","Sensitivity '$legadd'",true,"comparedcontours_th13dcp_datavssensitivity")'

mv comparedcontours_th13dcp_datavssensitivity.pdf $outfolder/comparedcontours_th13dcp_datavssensitivity.pdf
mv comparedcontours_th13dcp_datavssensitivity.png $outfolder/comparedcontours_th13dcp_datavssensitivity.png
mv comparedcontours_th13dcp_datavssensitivity.root $outfolder/comparedcontours_th13dcp_datavssensitivity.root

root -l -q -b 'utils/macros/makeOfficialPlotStyle.C++("'$outfolder'/comparedcontours_th13dcp_datavssensitivity.root","'$outfolder'/comparedcontours_th13dcp_datavssensitivity_official.root","T2K Run 1-7c preliminary","c1")'

#sort this bit
root -l -q -b 'utils/contours_app1D_comparison.C("'$datainputfile'","'$sensitivityinputfile'","Data '$legadd'","Sensitivity '$legadd'",0,true)'
mv contours_1D_dcp_compare.pdf $outfolder/contours_1D_dcp_datasensitivity_sincompare.pdf
mv contours_1D_dcp_compare.eps $outfolder/contours_1D_dcp_datasensitivity_sincompare.eps
mv contours_1D_dcp_compare.png $outfolder/contours_1D_dcp_datasensitivity_sincompare.png
mv contours_1D_dcp_compare.root $outfolder/contours_1D_dcp_datasensitivity_sincompare.root
root -l -q -b 'utils/macros/makeOfficialPlotStyle.C++("'$outfolder'/contours_1D_dcp_datasensitivity_sincompare.root","'$outfolder'/contours_1D_dcp_datasensitivity_sincompare_official.root","T2K Run 1-7c preliminary","c1")'

root -l -q -b 'utils/contours_app1D_comparison.C("'$datainputfile'","'$sensitivityinputfile'","Data '$legadd'","Sensitivity '$legadd'",0,false)'
mv contours_1D_dcp_compare.pdf $outfolder/contours_1D_dcp_datasensitivity_compare.pdf
mv contours_1D_dcp_compare.eps $outfolder/contours_1D_dcp_datasensitivity_compare.eps
mv contours_1D_dcp_compare.png $outfolder/contours_1D_dcp_datasensitivity_compare.png
mv contours_1D_dcp_compare.root $outfolder/contours_1D_dcp_datasensitivity_compare.root
root -l -q -b 'utils/macros/makeOfficialPlotStyle.C++("'$outfolder'/contours_1D_dcp_datasensitivity_compare.root","'$outfolder'/contours_1D_dcp_datasensitivity_compare_official.root","T2K Run 1-7c preliminary","c1")'