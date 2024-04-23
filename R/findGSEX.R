# findGSEX R Package - Estimate Genome Size
# Version 1.2
#
# Description:
# The findGSEX R package provides a method for estimating genome size by fitting
# k-mer frequencies in short reads with a normal distribution model.
#
# Author:
# Laiyi Fu & Hequan Sun
# Contact: laiyifu@xjtu.edu.cn or sunhequan@gmail.com
#
# Copyright (C) 2023-2024 Laiyi Fu & Hequan Sun (MPIPZ)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Installation:
# To install and use this package, source the provided R script:
# devtools::install_github("sperfu/findGSEX")
#
# Usage:
# Set options (optional):
# options(warn = -1)
#
# Define input parameters:
# path <- "histo_files"
# samples <- "O_21mer.histo"
# sizek <- 21
# exp_hom <- 200
# ploidy <- 4
# output_dir <- "outfiles"
# xlimit <- -1
# ylimit <- -1
# range_left <- exp_hom * 0.2
# range_right <- exp_hom * 0.2
#
# Call the findGSEX function with specified parameters:
# findGSEX(path, samples, sizek, exp_hom, ploidy, range_left, range_right, xlimit, ylimit, output_dir)
#
# For any questions, usage inquiries, or reporting potential bugs, please contact the author.
#
# License:
# This program is distributed under the terms of the GNU General Public License.
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <http://www.gnu.org/licenses/>.

##
## output:
##  *hist_hap_genome_size_est.pdf.
#


source('R/findGSE_v1.95_new.R')
library(scales)
library(png)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
## function for minimizing difference from original observed kmer-freq dr: mean, sd, xi, and scaling factor ##
start_time <- proc.time()
error_minimize<-function(tooptimize, x, end, xfit, xfit_left, xfit_right, d, min_valid_pos, itr)
{
  # x            : histogram replicated from d, which is formed from dr
  # xfit         : x-values to get the skew normal distribution
  # xfit-left    : a left-side  position to calculate initial mean and sd
  # xfit-right   : a right-side position to calculate initial mean and sd
  # d            : the observed kmer-freq that will be fitted
  # min_valid_pos: a left-side  position, from which the observed kmer-freq will be fitted
  # end          : a rigth-side position, till which the observed kmer-freq will be fitted
  sd_scale_factor    <- tooptimize[1]
  yfit_scale_factor  <- tooptimize[2]*100000/itr
  mean_scale_factor  <- tooptimize[3]
  xi_scale_factor    <- tooptimize[4]
  xfit               <- seq(min(x),max(x),length=end)
  meanfit            <- mean(x[x>=xfit_left & x<=xfit_right])*mean_scale_factor
  sdfit              <- sd(x[  x>=xfit_left & x<=xfit_right])*sd_scale_factor
  xifit              <- 1*xi_scale_factor
  # cat(paste("   Info: initialized mean as ", meanfit, "\n", sep=""))
  # cat(paste("   Info: initialized sd   as ", sdfit, "\n", sep=""))
  #cat(paste("   Info: initialized skew as ", xifit, "\n", sep=""))
  #yfit               <- dsnorm(xfit,mean=meanfit,sd=sdfit, xi=xifit)*yfit_scale_factor*length(x)
  yfit               <- dnorm(xfit,mean=meanfit,sd=sdfit)*yfit_scale_factor*length(x)
  diff0              <- sqrt(sum((d[min_valid_pos:end, 2] - yfit[min_valid_pos:end])^2))
  #diff0              <- sqrt(sum((d[xfit_left:end, 2] - yfit[xfit_left:end])^2))
  # cat(paste("   111 Info: initialized min_valid_pos   as ", min_valid_pos, "\n", sep=""))
  # cat(paste("   111 Info: initialized end as ", end, "\n", sep=""))
  # cat(paste("   111 Info: initialized x_fit_left   as ", xfit_left, "\n", sep=""))
  # cat(paste("   111 Info: initialized x_fit_right as ", xfit_right, "\n", sep=""))
  # diff0              <- sqrt(sum((d[xfit_left:xfit_right, 2] - yfit[xfit_left:xfit_right])^2))
  #
  return(diff0)
}
# function for tunning final fitting if heterozgyous genomes
error_minimize2<-function(tooptimize, h_het, h_hom, h_target)
{
  # h_het   : raw fitting for the heterozygous region
  # h_hom   : raw fitting for the homozygous   region
  # h_target: before valley: fitted cnt, after valley: raw cnt
  het_delta <- 1-tooptimize[1]
  hom_delta <- 1 # reserved parameter: might be applied in future
  h_merge   <- het_delta*h_het + hom_delta*h_hom
  diff2     <- sum(abs(h_merge - h_target))
  #
  return(diff2)
}
## function for minimizing difference from remainder kmer-freq dr: mean, sd, and scaling factor ##
error_minimize3<-function(tooptimize, x, end, xfit_left, xfit_right, d, min_valid_pos, itr)
{
  # x            : histogram
  # xfit         : x-values to get the normal distribution
  # xfit-left    : a left-side  position to calculate initial mean and sd
  # xfit-right   : a right-side position to calculate initial mean and sd
  # d            : the observed kmer-freq that will be fitted
  # min_valid_pos: a left-side  position, from which the observed kmer-freq will be fitted
  # end          : a rigth-side position, till which the observed kmer-freq will be fitted
  sd_scale_factor    <- tooptimize[1]
  yfit_scale_factor  <- tooptimize[2]*100000/itr
  mean_scale_factor  <- tooptimize[3]
  xfit               <- seq(min(x),max(x),length=end)
  meanfit            <- mean(x[x>=xfit_left & x<=xfit_right])*mean_scale_factor
  sdfit              <- sd(x[  x>=xfit_left & x<=xfit_right])*sd_scale_factor
  yfit               <- dnorm(xfit,mean=meanfit,sd=sdfit)*yfit_scale_factor*length(x)
  diff0              <- sqrt(sum((d[min_valid_pos:end, 2] - yfit[min_valid_pos:end])^2))
  # cat(paste("   111 Info: initialized min_valid_pos   as ", min_valid_pos, "\n", sep=""))
  # cat(paste("   111 Info: initialized end as ", end, "\n", sep=""))
  # cat(paste("   111 Info: initialized x_fit_left   as ", xfit_left, "\n", sep=""))
  # cat(paste("   111 Info: initialized x_fit_right as ", xfit_right, "\n", sep=""))
  # diff0              <- sqrt(sum((d[xfit_left:xfit_right, 2] - yfit[xfit_left:xfit_right])^2))
  #diff0              <- sqrt(sum((d[xfit_left:end, 2] - yfit[xfit_left:end])^2))
  #
  return(diff0)
}
# recover count 0: in initial count, making kmer freq consecutive.
initial_count_recover <- function(d0)
{
  # dr: the initial kmer count from softwares like jellyfish
  dr <- cbind(1:d0[length(d0[,1]), 1], rep(0, d0[length(d0[,1]), 1]))
  i  <- 0
  while( i < length(d0[,1]))
  {
    i               <- i + 1
    dr[d0[i, 1], 2] <- d0[i, 2]
  }
  return (dr)
}
# modify kmer freq before fitting
kmer_count_modify <- function(start, end, left_right, histx)
{
  # start or end does not include the peak position
  # modify rigth part according to left
  if(left_right == 0)
  {
    for (i in start:end)
    {
      histx[2*end+2-i, 2]   <- histx[i, 2]
    }
  }else
  {
    for (i in start:end)
    {
      histx[2*start-2-i, 2] <- histx[i, 2]
    }
  }
  return (histx)
}
#
findGSE_sp <- function(histo="", sizek=0, outdir="", exp_hom=0, species="",ploidy_ind=2,avg_cov = 0,left_fit_ratio = 0.835)
{
  # initial values
  if(missing(histo))
  {
    stop(paste("Cannot find histo file: ", histo, ". Program exited.\n",sep=""))
  }
  if(missing(sizek))
  {
    stop(paste("Cannot find sizek info: ", sizek, ". Program exited.\n",sep=""))
  }
  # defaults
  cat('ploidy index start: ',ploidy_ind,'\n')
  if(missing(exp_hom))  exp_hom <- 0
  if(missing(outdir))   outdir  <- getwd()
  ######################################## libararies required ###############################################
  #cat("\n special version: \n")
  # find all peaks in a time series
  suppressWarnings(suppressMessages(library("pracma")))
  # provide skew normal distrbution fitting (dsnorm(...))
  suppressWarnings(suppressMessages(library("fGarch")))
  ## version id
  vers <- 'v1.94.'
  ##
  cat("\nfindGSE initialized...\n")
  ## running time
  start_t <- Sys.time()
  ## check args
  # expert-tunning het fitting (later it would be automatically adjusted)
  het_fitting_delta <- 0.11
  histof  <- F
  kpro    <- F
  if(histo != "") histof <- T
  if(sizek >  0)  kpro   <- T
  if(histof == F || kpro == F)
  {
    stop("Option -histo or -sizek was provided incorrectly. Program exited.\n")
  }
  ## default args
  countingfile    <- histo
  if(file.exists(countingfile)==F)
  {
    stop(paste("Cannot find histo file: ", countingfile, ". Program exited.\n",sep=""))
  }
  cat(paste("   Info: histo file provided as ", countingfile, "\n", sep=""))
  #
  targetsizek     <- sizek
  if(targetsizek  < 15)
  {
    stop(paste("Unexpected size of k: ", targetsizek, ". Size k>=15 required. Program exited.\n",sep=""))
  }
  cat(paste("   Info: size k set as ", targetsizek, "\n", sep=""))
  #
  path_set        <- F
  path            <- outdir
  if(path==".")
  {
    path <- getwd()
  }else
    if(path == "")
    {
      path         <- dirname(normalizePath(countingfile))
    }
  path_set       <- T
  cat(paste("   Info: output folder set as ", path, "\n", sep=""))
  #
  het_observed    <- F
  # expected maximum value for the hom k-mer peak! set as fpeak~2*fpeak
  exp_hom_cov_max <- exp_hom
  cat(paste("   Info: expected coverage of homozygous k-mers set as ", exp_hom_cov_max, "\n", sep=""))
  exp_hom_cov_max_save <- exp_hom_cov_max
  # derived  maximum value for the het k-mer peak!
  der_het_cov_max <- 0
  ## if output folder does not exist, create it; if existing, give a warning.
  dir.create(file.path(path), showWarnings = F)
  if(exp_hom_cov_max > 0)
  {
    het_observed      <- T
    # derived  maximum value for the het k-mer peak.
    der_het_cov_max   <- round(ploidy_ind/(ploidy_ind+1)*exp_hom_cov_max);
    cat(paste("   Info: het observed set as true ==> heterozygous fitting asked. \n", sep=""))
    if(het_fitting_delta != 0.11)
    {
      cat(paste("   Info: value for tuning het fitting set as ", het_fitting_delta, "\n", sep=""))
    }
  }else
  {
    het_observed      <- F
    only_one_het_peak <- F
    only_one_hom_peak <- T
    cat(paste("   Info: het observed set as false ==> heterozygous fitting not asked. \n", sep=""))
  }
  ##
  prefix <- ''
  samples          <- basename(countingfile)
  for (sample in samples)
  {
    cat(paste("\nIterative fitting process for sample ", sample, " started... \n", sep=""))
    # all accessions -- where i can change and modify
    # initialize variables where outputs are recorded
    size_all       <-  0   # Gsize with               error-k-mer-freq
    size_exl       <-  0   # Gsize without            error-k-mer-freq
    size_fit       <-  0   # Gsize with              fitted-k-mer-freq
    size_cat       <-  0   # Gsize with fitted-cat-observed-k-mer-freq
    size_cor       <-  0   # Gsize corrected with size_cat and skewness; size_cor=size_cat/skewness (no)
    size_cor2      <-  0   # Gsize corrected with fitted-cat-observed-freq and upperbound freq 'e' (yes)
    i              <-  0   # the ith sample/accession
    kmer_peak_cnt  <-  0   # k-mer count at the peak position of k-mer counting
    kmer_peak_pos  <-  0   # genome-wide k-mer coverage
    kmer_peak_pos2 <-  0   # genome-wide k-mer coverage according to fitted curve
    labels         <- NA   # labels according input below
    namebank       <- NA
    border_pos     <-  1
    min_valid_pos  <-  1
    first_peak_pos <-  1
    first_mean     <-  0
    first_sd       <-  1
    first_mean_raw <-  0
    # initialize parameters for several iterations of fitting
    mymeanfit      <-  0
    mysdfit        <-  1
    myxifit        <-  1
    myscalefit     <-  1
    # special fit at itr 1
    norm_fit_itr1  <-  F
    # open pdf to write visualization of fitting
    pdf(paste(path, '/',prefix, vers, 'est.', sample, '.sizek', targetsizek, '.curvefitted.pdf',
              sep=""))
    first_fit  <- TRUE
    for (sizek in targetsizek[1]) {
      normal   <- 0
      if(file.exists(countingfile)==TRUE)
      {
        dhet   <- read.table(countingfile)
        ## recover positions in dr whose initial counts are 0
        dr     <- initial_count_recover(dhet[1:length(dhet[,1]),])
        rm("dhet")
        dhet   <- dr
        # if heterozygosity oberved in k-mer freq,
        #    fit het and remove fitted values from subsequent hom-fitting
        selected    <- min(2000, length(dhet[,1]))
        if(selected < exp_hom_cov_max)
        {
          selected  <- min(round(1.5*exp_hom_cov_max), length(dhet[,1]))
        }
        main_peak_is_hom <- T
        if(het_observed)
        {
          # initialize count as error for fitting
          error    <- dr
          cat(paste("    Size ", sizek , " fitting for het k-mers \n", sep=""));
          # find hom and het peaks (if any).
          peaks    <- findpeaks(error[2:exp_hom_cov_max, 2])
          #peaks    <- filter_peaks(tmp_peaks, dhet)
          if(length(peaks) == 0)
          {
            plot(dhet[, 1], dhet[, 2],
                 col="gray",
                 xlab="K-mer freq", ylab="Number of k-mers",
                 xlim=c(1, 4*(which.max(dhet[5:selected, 2])+4)), ylim=c(1, 1.5*max(dhet[5:selected, 2])))
            text(1.5*(which.max(dhet[5:selected, 2])+4), 1.1*max(dhet[5:selected, 2]),
                 paste("Main peak (after freq=5) at k-mer freq: ",
                       which.max(dhet[5:selected, 2])+4,
                       sep=""))
            dev.off();
            stop("No k-mer freq peak was found with the cutoff given by -exp_hom.
                  You may need to increase the cutoff!")
            quit("no")
          }
          hom_peak_pos <- 0
          het_peak_pos <- 0
          cat('enter ploidy equal to ',ploidy_ind)
	  if(ploidy_ind == 1){
          hom_peak_pos <- peaks[which.max(peaks[,1]), 2]+1   # hom peak
          peaks[which.max(peaks[,1]), 1] <- 0
          het_peak_pos <- peaks[which.max(peaks[,1]), 2]+1   # het peak
          peaks[which.max(peaks[,1]), 1] <- 0
	  }else if(ploidy_ind == 10){

	  het_peak_pos <- 462 #ploidy_ind * avg_cov
	  hom_peak_pos <- 615 #(ploidy_ind + 1) * avg_cov
	  }
	  else{
	  het_peak_pos <- ploidy_ind * avg_cov
	  hom_peak_pos <- (ploidy_ind + 1) * avg_cov
	  }
	  cat('hom_peak_pos is: ',hom_peak_pos,'\n')
	  cat('het_peak_pos is: ',het_peak_pos,'\n')
          # if two peaks found, the larger one is hom; the other is het
          if(het_peak_pos > hom_peak_pos)
          {
            tmp               <- het_peak_pos
            het_peak_pos      <- hom_peak_pos
            hom_peak_pos      <- tmp
            main_peak_is_hom  <- F
            cat("    Warning: two peaks found in k-mer freq with given -exp_hom,
                  peak with lower height as hom-peak!\n")
          }
          if(het_peak_pos>0 & hom_peak_pos>0)
          {
            # in case two peask are near each other (e.g. it may happen with k-mer cov >100x)
            ratio        <- hom_peak_pos/het_peak_pos;
            if (round(ratio)<1.5) # nearly equal
            {
              hom_peak_pos <- max(hom_peak_pos, het_peak_pos)
              het_peak_pos <- min(hom_peak_pos, het_peak_pos)
            }
          }
          # if one peak  found, need to check if it is hom or het accordingly
          only_one_hom_peak <- F
          only_one_het_peak <- F
          if(het_peak_pos == hom_peak_pos)
          {
            unkownpeak    <- het_peak_pos
            if(unkownpeak >  der_het_cov_max)
            {
              het_peak_pos      <- round(ploidy_ind/(ploidy_ind+1)*hom_peak_pos)
              cat("    Warning: only one peak in k-mer freq with given -exp_hom,
                  determined as hom-peak!\n")
              cat("            -exp_hom needs to be increased
                  if you think it is a het peak.\n")
              only_one_hom_peak <- T
            }else
              if(unkownpeak  < der_het_cov_max)
              {
                hom_peak_pos      <- (ploidy_ind+1)/ploidy_ind*het_peak_pos
                cat("    Warning: only one peak in k-mer freq with given -exp_hom,
                    determined as het-peak!\n")
                cat("          -exp_hom needs to be decreased
                    if you think it is a hom peak.\n")
                main_peak_is_hom  <- F
                only_one_het_peak <- T
              }else
              {
                plot(dhet[, 1], dhet[, 2],
                     col="gray",
                     xlab="K-mer freq", ylab="Number of k-mers",
                     xlim=c(1, selected), ylim=c(1, 1.5*max(dhet[3:selected, 2])))
                text(100, max(dhet[3:selected, 2]),
                     paste("Main peak at k-mer freq: ",
                           which.max(dhet[3:selected, 2])+2,
                           sep=""))
                dev.off()
                stop("Unexpected peaks found -
                     please check the raw k-mer curve provided in pdf,
                     and reset the option -exp_hom with a proper value x,
                     which must satisfy 2*hom_peak>x>hom_peak !\n\n");
              }
          }
          cat('    Info: het_peak_pos for het fitting: ', het_peak_pos, '\n')
          cat('    Info: hom_peak_pos for hom fitting: ', hom_peak_pos, '\n')
          het_peak_pos_save <- het_peak_pos
          hom_peak_pos_save <- hom_peak_pos
          het_min_valid_pos <- which.min(error[1:het_peak_pos, 2])
          if(ploidy_ind == 1){
          # fit het distribution
          het_xfit_left     <- max(min(het_min_valid_pos, het_peak_pos-10), 3)        # update on 20200129
          if(het_peak_pos   > 80)
          {
            het_xfit_left  <- floor(het_xfit_left+het_peak_pos*0.1473684)
          }
          het_xfit_right   <- het_peak_pos + round(0.5*(hom_peak_pos-het_peak_pos))   # update on 20200129
          if(het_peak_pos   > 80)
          {
            het_xfit_right <- floor(het_xfit_right-het_peak_pos*0.1473684)
          }
	  }else if(ploidy_ind == 10){

          het_xfit_left <- 350 #ploidy_ind * avg_cov
          het_xfit_right <- 650 #(ploidy_ind + 1) * avg_cov
          }
          else{
	    het_xfit_left <- het_peak_pos * left_fit_ratio  # 0.85  0.83
	    het_xfit_right <- het_peak_pos * 1.15  # 1.15  1.13
	  }
          if(het_xfit_right > het_peak_pos-1 - het_xfit_left + het_peak_pos-1)
          {
            het_xfit_right  <- het_peak_pos-1 - het_xfit_left + het_peak_pos-1
          }
          histx2    <- kmer_count_modify(het_xfit_left, het_peak_pos-1, 0, error)
          x         <- rep(1:het_xfit_right,ceiling(abs(histx2[1:het_xfit_right,2])/(100000/1)))
          xfit      <- seq(min(x),max(x),length=het_xfit_right)
          # optimizing parameters
          cat('    Info: het_xfit_left  for het fitting: ', het_xfit_left, '\n')
          cat('    Info: het_xfit_right for het fitting: ', het_xfit_right, '\n')
          het_xfit_right_save <- het_xfit_right
          v<-optim(par=c(1, 1, 1, 0.8),
                   fn=error_minimize, x=x, end=het_xfit_right, xfit=xfit,
                   xfit_left=het_xfit_left, xfit_right=het_xfit_right,
                   d=histx2, min_valid_pos=het_min_valid_pos, itr=1)
          #
          # scale and correct yfit with optimized values in v$par
          xfit2     <- seq(0.01,selected,0.01)
          meanfit   <- mean(x[x>=het_xfit_left & x<=het_xfit_right])*v$par[3]
          if(meanfit < 0)
          {
            cat(paste("    Warning: data does not follow
                      assumed distribution anymore; fitting stopped.\n", sep=""))
          }
          sdfit     <- sd(x[  x>=het_xfit_left & x<=het_xfit_right])*v$par[1] # raw sd
          #
          xifit     <- 1*v$par[4];
          #yfit0     <- dsnorm(xfit2,mean=meanfit,sd=sdfit, xi=xifit)*(100000/1)*v$par[2]*length(x)

          yfit0     <- dnorm(xfit2,mean=meanfit,sd=sdfit)*(100000/1)*v$par[2]*length(x)
          #
          ## scale fitting according to observation
          hetfit    <- yfit0[seq(100, length(yfit0), 100)]
          hetfit    <- hetfit*dhet[het_peak_pos,2]/max(hetfit[1:selected])
          ## end scaling ##
          #
          # set d0 for further hom-fitting
          d0        <- dr
          if(only_one_hom_peak)
          {
            d0[1:selected , 2] <- d0[1:selected , 2] - (1-het_fitting_delta)*hetfit[1:selected]
          }else
          {
            d0[1:selected , 2] <- d0[1:selected , 2] - hetfit[1:selected]
          }
        }else
        {
          d0                 <- dr
          hetfit             <- 0
          hetfit[1:selected] <- 0
          het_min_valid_pos  <- 100000
          only_one_hom_peak  <- T
        }
        i    <- i + 1
        if(i > length(targetsizek)) { break }
        rm("dr")
        # dr can be different from raw k-mer counting due to het fit
        dr   <- d0
        #
        # initialize count as error for fitting: initial error vector
        error        <- dr
        # set the right-bound for fitting
        maxright4fit <- dr[length(dr[,1])-1, 1]
        maxright4fit <- min(10000,length(dr[,1]))
        # set the number of total iterations corresponding
        #     to raw count given by softwares such as jellyfish
        totalitr     <- 2
        myitr        <-  0
        myyfit       <-  0
        homfit       <-  0
        #
        first_fit    <- TRUE
        itr = 1
        while(itr <= totalitr)
        {
          if(first_fit)
          {
            cat(paste("    Size ", sizek , " at itr ", itr, "\n", sep=""))
            if(itr>=2 && 2*mymeanfit[itr-1] > maxright4fit) { myitr <- itr - 1; break;}
            if(itr == 1) # original data
            {
              # 1st-col    : height,
              # 2nd-col    : index of maximum,
              # 3rd/4th-col: indices that peak begins/ends
              if(het_observed)
              {
                # for  border_pos
                peaks           <- findpeaks(error[3:selected, 2])
                peaks[,2:4]     <- peaks[,2:4]+2
                border_pos      <- peaks[1, 3]          # caution: peaks[peakindex, 3]
                # for others: caution!
                if(main_peak_is_hom==F)
                {
                  peaks           <- findpeaks(error[round(0.5*(het_peak_pos+hom_peak_pos)):selected, 2])
                  peaks[,1:4]     <- peaks[,1:4]+round(0.5*(het_peak_pos+hom_peak_pos))-1
                }else
                {
                  peaks           <- findpeaks(error[het_peak_pos:selected, 2])
                  peaks[,1:4]     <- peaks[,1:4]+het_peak_pos-1
                }
              }
              else
              {
                peaks           <- findpeaks(error[3:selected, 2])
                peaks[,2:4]     <- peaks[,2:4]+2
                ## new on 2017-09-16
                peaks           <- filter_peaks(peaks, error)
                ## end 2017-09-16
                border_pos      <- peaks[1, 3]          # caution: peaks[peakindex, 3]
              }
              peakindex         <- which.max(peaks[,1]) # caution: if double peaks
              min_valid_pos     <- peaks[1, 3]          # caution: peaks[peakindex, 3]
              original_peak_pos <- peaks[peakindex, 2]
              original_peal_val <- peaks[peakindex, 1]
              if(het_observed && only_one_het_peak)
              {
                original_peak_pos <- (ploidy_ind+1)/ploidy_ind*het_peak_pos
                original_peal_val <- peaks[original_peak_pos, 1]
                min_valid_pos     <- het_peak_pos
              }
              datacheck         <- min(300, maxright4fit)
              dci               <- 0
              cat('    Info: min_valid_pos: ', min_valid_pos,'\n')
              min_valid_pos_save = min_valid_pos
              while(min_valid_pos == original_peak_pos)
              {
                min_valid_pos     <- min_valid_pos + 1
                original_peak_pos <- which.max(error[min_valid_pos:maxright4fit,2])+(min_valid_pos-1)
                original_peal_val <- max(error[min_valid_pos:maxright4fit,2])
                dci               <- dci + 1
                if(dci > datacheck) { normal=1;break; }
              }
              if(dci   > datacheck)
              {
                cat(paste("Error: weird initial count-data of sample ",
                          sample, ". Please check!!!\n",
                          sep=""))
                normal <- 1
                break
              }
              cat('    Info: signal error border: ', border_pos,'\n')
              first_peak_pos  <- original_peak_pos
              # copy
              d               <- error
              # complement the left side of the peak position
              #    according to observations on the right side of the peak
              for (li in 1:(original_peak_pos-1))
              {
                d[li, 2] <- d[min(2*original_peak_pos-li, length(dr[,1])),2]
              }
              end <- min(original_peak_pos+round((original_peak_pos-min_valid_pos)*5/12),length(dr[,1]))
              if(end - original_peak_pos <= 3 && original_peak_pos>=10)
              {
                end      <- end + 5
              }
            }
            else # error reminder
            {
              # caution: valid kmer count is at least 3
              min_valid_pos     <- round(mymeanfit[itr-1])+3
              original_peak_pos <- which.max(error[min_valid_pos:maxright4fit, 2])+(min_valid_pos-1)
              original_peal_val <- max(error[min_valid_pos:maxright4fit, 2])
              datacheck         <- min(300, maxright4fit)
              dci               <- 0
              while(min_valid_pos == original_peak_pos)
              {
                min_valid_pos     <- min_valid_pos + 1
                original_peak_pos <- which.max(error[min_valid_pos:maxright4fit,2])+(min_valid_pos-1)
                original_peal_val <- max(error[min_valid_pos:maxright4fit,2])
                dci               <- dci + 1
                if(dci > datacheck) { break; }
              }
              if(dci   > datacheck)
              {
                cat(paste("    Warning: weired initial count-data of sample ",
                          sample, ". Please check!!!\n", sep=""))
                normal   <- 0
                xlimmax  <- 100
                break
              }
              # copy
              d          <- abs(error)
              # complement the left side of the peak position according to observations on the right side of the peak
              thisnormal <- 0
              for (li in (original_peak_pos-round(mymeanfit[1])):(original_peak_pos-1)) {
                if(2*original_peak_pos-li > maxright4fit){ thisnormal<-1; break}
                d[li, 2] <- d[2*original_peak_pos-li,2]
              }
              if(thisnormal==1){normal=0; break;}
              end        <- min(original_peak_pos+round(mymeanfit[itr-1]), maxright4fit)
              if(end - original_peak_pos <= 3 && original_peak_pos>=10)
              {
                end      <- end + 5
              }
            }
          }
          else
          {
            # due to the large diff in the first fitting and original curve, increase end for fitting.
            end           <- min(2*original_peak_pos-1, maxright4fit)
            min_valid_pos <- min_valid_pos + 1 ##### 2016-09-16 new
            d             <- error
            first_fit     <- TRUE
          }
          # prepare the histogram data with replication
          forrep     <- ceiling(abs(d[1:end,2])/(100000/itr))
          anyNA      <- which(is.na(forrep))
          if (length(anyNA) >= 1)
          {
            d[anyNA, 2] <- dr[anyNA, 2]; # caution: missing values could happen here!!!
            cat(paste("    Warning: ", length(anyNA),
                      " data-points of d has been replaced as 1 to properly prepare hisgram data x
                      for function optim at itr ", itr, "\n", sep=""))
          }
          x          <- rep(1:end,ceiling(abs(d[1:end,2])/(100000/itr)))
          xfit       <- seq(min(x),max(x),length=end)
          # normal fit
          matchset   <- '0'
          if(!is.na(match(sizek,matchset)))
          {
            xfit_left  <- min_valid_pos
            xfit_right <- 2*original_peak_pos-xfit_left
            if(first_peak_pos>=200)
            {
              #xfit_left  <- original_peak_pos - 50
              #xfit_right <- original_peak_pos + 50
              xfit_left  <- original_peak_pos * left_fit_ratio #- 50 0.85
              xfit_right <- original_peak_pos * 1.13 #+ 50 1.15
            }
            v<-optim(par=c(1, 1, 1, 1),
                     fn=error_minimize, x=x, end=end, xfit=xfit,
                     xfit_left=xfit_left, xfit_right=xfit_right,
                     d=error,
                     min_valid_pos=min_valid_pos, itr=itr)
          }
          else
          {
            xfit_left  <- min_valid_pos
            xfit_right <- 2*original_peak_pos-1
            if(first_peak_pos>=100 && itr<2)
            {
              xfit_left  <- original_peak_pos - 20
              xfit_right <- original_peak_pos + 20
            }else
              if(first_peak_pos>=30  && itr<2)
              {
                xfit_left  <- original_peak_pos - 10
                xfit_right <- original_peak_pos + 10
              }
            v<-optim(par=c(1, 1, 1, 0.8),
                     fn=error_minimize, x=x, end=end, xfit=xfit,
                     xfit_left=xfit_left, xfit_right=xfit_right,
                     d=error,
                     min_valid_pos=min_valid_pos, itr=itr)
          }
          cat('    Info: hom_xfit_left  for hom fitting at itr ', itr,  ': ', xfit_left, '\n')
          cat('    Info: hom_xfit_right for hom fitting at itr ', itr,  ': ', xfit_right, '\n')
          # if(itr == 1){
          #   xfit_left
          # }
          #
          # scale and correct yfit with optimized values in v$par
          xfit2     <- seq(0.01,maxright4fit,0.01)
          meanfit   <- mean(x[x>=xfit_left & x<=xfit_right])*v$par[3]
          if(meanfit < 0)
          {
            cat(paste("    Warning: data does not follow assumed distribution anymore at itr ",
                      itr, ", fitting stopped.\n", sep=""))
            normal  <- 0
            if(myyfit==0 && itr==1)
            {
              v<-optim(par=c(1, 1, 1),
                       fn=error_minimize3, x=x, end=end,
                       xfit_left=xfit_left, xfit_right=xfit_right,
                       d=error,
                       min_valid_pos=min_valid_pos, itr=itr)
              xfit2    <- seq(0.01,maxright4fit,0.01)
              yfit0    <- dnorm(xfit2,
                                mean=meanfit,
                                sd=sdfit)*(100000/itr)*v$par[2]*length(x)
              myyfit   <- myyfit + yfit0
              norm_fit_itr1 <- T
              if(het_observed && itr == 1)
              {
                homfit  <- yfit0
              }
            }
            break
          } # stop if the data does not follow normal distribution anymore
          sdfit     <- sd(x[  x>=xfit_left & x<=xfit_right])*v$par[1]
          xifit     <- 1*v$par[4]
          #yfit0     <- dsnorm(xfit2,mean=meanfit,sd=sdfit, xi=xifit)*(100000/itr)*v$par[2]*length(x)
	        yfit0     <- dnorm(xfit2,mean=meanfit,sd=sdfit)*(100000/itr)*v$par[2]*length(x)
          yfit      <- yfit0[seq(100, length(yfit0), 100)]
          ## repeat first fitting if fitting is not satisfactory
          if(itr==1 && end<min(2*original_peak_pos-1, maxright4fit))
          {
            fitting_check <- error[1:maxright4fit, 2] - yfit
            n_less_than_0 <- 0
            for(fci in c((original_peak_pos+3):min(round(2*original_peak_pos), length(dr[,1]))))
            {
              if(fitting_check[fci] < 0)
              {
                n_less_than_0 <- n_less_than_0 + 1
              }
            }
            if(n_less_than_0 >= 3)
            {
              first_fit <- FALSE
              cat(paste("       Fitting has to be repeated: sizek ", sizek,
                        " at itr ", itr, "\n",
                        sep=""))
              i
              next
            }
          }
          # plot
          if(i == 1)
          {
            #plot(error[1:maxright4fit, 1], abs(error[1:maxright4fit, 2]),
            #      xlim=c(0, min(maxright4fit, 4*meanfit)), ylim=c(0, 1.5*max(error[xfit_left:maxright4fit, 2])),
            #      type='b', col="gray",
            #      xlab="K-mer Coverage", ylab="Genomewide k-mer Frequency")
            # title(main=paste('Example: fitting at itr ', itr, 'of sample ', sample, 'k=', sizek, ')'),
            #     cex.main=0.8)
            # lines(xfit2[100:length(xfit2)], yfit0[100:length(xfit2)], col="blue", lwd=2)
          }
          #### begin
          if(max(yfit) > 1.2*max(error[min_valid_pos:maxright4fit, 2]))
          {
            cat(paste("    Note on hom fitting: fitting stopped at iter ",
                      itr, ", expected: ", totalitr, "\n",sep=""))
            if(myyfit==0 && itr==1)
            {
              v<-optim(par=c(1, 1, 1),
                       fn=error_minimize3, x=x, end=end,
                       xfit_left=xfit_left, xfit_right=xfit_right,
                       d=error,
                       min_valid_pos=min_valid_pos, itr=itr)
              xfit2    <- seq(0.01,maxright4fit,0.01)
              yfit0    <- dnorm(xfit2,
                                mean=meanfit,
                                sd=sdfit)*(100000/itr)*v$par[2]*length(x)
              if(sum(yfit0) > 0)
              {
                myyfit   <- myyfit + yfit0
                norm_fit_itr1 <- T
                if(het_observed && itr == 1)
                {
                  homfit  <- yfit0
                }
              }else # abnormal case: caution!: no fitting performed on hom fitting
              {
                yfit0    <- rep(0, length(xfit2))
                myyfit   <- myyfit + yfit0
                norm_fit_itr1 <- T
                if(het_observed && itr == 1)
                {
                  homfit  <- yfit0
                }
              }
            }
            break;
          }
          #### end
          # estimate error
          error[1:maxright4fit, 2] <- abs(error[1:maxright4fit, 2] - yfit)
          # get more accurate fitted values
          xfit2    <- seq(0.01,maxright4fit,0.01);
          #yfit0    <- dsnorm(xfit2,mean=meanfit,sd=sdfit, xi=xifit)*(100000/itr)*v$par[2]*length(x)


          yfit0    <- dnorm(xfit2,mean=meanfit,sd=sdfit)*(100000/itr)*v$par[2]*length(x)
          myyfit   <- myyfit + yfit0
          # record fitted parameters
          mymeanfit[itr]  <- meanfit
          mysdfit[itr]    <- sdfit
          myxifit[itr]    <- xifit
          myscalefit[itr] <- (100000/itr)*v$par[2]*length(x)
          if(itr == 1)
          {
            first_mean     <- meanfit
            first_sd       <- sdfit
          }
          #
          ## het est: new feature
          if(het_observed && itr == 1)
          {
            homfit  <- yfit0
          }
          ## end of het rate est: new feature
          #
          itr <- itr + 1
        }
        cat("Iterative fitting done.\n\n")
        if(normal == 0)
        {
          # original count
          xlimmax <- min(3*first_peak_pos, selected)
          plot(dhet[, 1], abs(dhet[, 2]),
               xlim=c(0,xlimmax), ylim=c(0, 1.5*max(dhet[border_pos:xlimmax, 2])),
               type='b', col="gray",
               xlab="K-mer freq", ylab="Number of k-mers")
          title(main=paste('Sample ', sample, ' k=', sizek), cex.main=0.8) #
          # fitted count
          if(sum(hetfit[1:selected] > 0)) # plus het
          {
            # fitted peak
            fitted_peak_value <- which.max(myyfit)/100
            if(norm_fit_itr1 == T)
            {
              fitted_peak_value <- hom_peak_pos
            }
            fittedyvalues <- myyfit[seq(100, length(myyfit), 100)]

            if(only_one_hom_peak)
            {
              # begin find optinum het_fitting_delta 2016-12-20
              h_target <- c(hetfit[1:min(border_pos, het_min_valid_pos)],
                            as.vector(dhet[min(border_pos, het_min_valid_pos)+1:(length(dr[,2])-min(border_pos, het_min_valid_pos)),2]))
              v2 <- optim(par=c(0.11, 1),
                          fn=error_minimize2,
                          h_het=hetfit[1:first_peak_pos],
                          h_hom=fittedyvalues[1:first_peak_pos],
                          h_target=h_target[1:first_peak_pos])
              ## update delta
              het_fitting_delta <- v2$par[1]
              cat(paste("    Info: het_fitting_delta optimized as ", het_fitting_delta, "\n\n", sep=""))
              #cat(paste("    Info: hom_fitting_delta optimized as ", v2$par[2],         "\n\n", sep=""))
              # end   find optinum het_fitting_delta 2016-12-20
            }

            if(only_one_hom_peak)
            {
              fittedyvalues[1:selected] <- fittedyvalues[1:selected] + (1-het_fitting_delta)*hetfit[1:selected]
            }else
            {
              fittedyvalues[1:selected] <- fittedyvalues[1:selected] + hetfit[1:selected]
            }
            lines(1:length(fittedyvalues), fittedyvalues, col="blue", lwd=5)
            # het
            if(only_one_hom_peak)
            {
              lines(1:selected, (1-het_fitting_delta)*hetfit[1:selected], col="cyan", lwd=2, lty="dashed")
            }
            else
            {
              lines(1:selected,  hetfit[1:selected], col="cyan", lwd=2, lty="dashed")
              options(scipen=999)
              het_fit <- cbind(1:selected,  round(hetfit[1:selected]))
              write.table(het_fit, file = paste(path, '/',prefix, vers, 'est.', sample, ".genome.size.estimated.k", min(targetsizek), 'to', max(targetsizek),".fitted_hetfit_count.txt", sep=""), sep = " ", quote=F, row.names = F, col.names = F)
	      cat(paste0('writing hetfit_count file done!!!\n het_fit length: ',dim(het_fit)[1]," lines.\n"))
              options(scipen=0)
            }
          }
          else
          {
            lines(xfit2[100:length(xfit2)], myyfit[100:length(xfit2)], col="blue", lwd=5)
            # fitted peak
            fitted_peak_value <- which.max(myyfit)/100
            fittedyvalues     <- myyfit[seq(100, length(myyfit), 100)]
          }
          # error rate of missing fitting: from raw k-mer counting dhet which includes het k-mers
          sum(fittedyvalues[border_pos:maxright4fit]-dhet[border_pos:maxright4fit,2])/sum(dhet[border_pos:maxright4fit,2])
          # fitted catting original with border_pos <- min_valid_pos; from raw k-mer counting dhet which includes het k-mers
          yfit2 <- c(fittedyvalues[1:min(border_pos, het_min_valid_pos)],
                     as.vector(dhet[min(border_pos, het_min_valid_pos)+1:(length(dr[,2])-min(border_pos, het_min_valid_pos)),2]))
          #
          lines(dhet[, 1], yfit2, col="red", lwd=1.5)
          #title(main=paste('Sample ', sample, ' k=', sizek, ')'), cex.main=0.8)
          #
          # calculate genome size (Mb) according to original data
          # with all kmers
          genome_size                 <- round(sum(dr[,1]*dhet[,2]/first_peak_pos))
          # with valid kmers: caution: 1.peak-value 2. valid kmer with cov >= min_valid_pos (or border_pos)
          genome_size_filtering_error <- round(sum(dr[border_pos:length(dhet[,1]),1]*dhet[border_pos:length(dr[,1]),2]/first_peak_pos))
          ##
          if(het_observed)
          {
            genome_size                 <- round(sum(dr[,1]*dhet[,2]/hom_peak_pos))
            genome_size_filtering_error <- round(sum(dr[border_pos:length(dhet[,1]),1]*dhet[border_pos:length(dr[,1]),2]/hom_peak_pos))
          }
          # genomic kmers
          # cat(paste("        err-exl genomic ", sizek, "mers: ", sum(dhet[border_pos:length(dhet[,2]),2]), " Mio\n", sep=""));
          # fitted data with fitted peak
          genome_size_fitted <- round(sum(c(1:length(fittedyvalues))*fittedyvalues[1:length(fittedyvalues)]/fitted_peak_value))
          # genomic kmers
          # cat(paste("        fit     genomic ", sizek, "mers: ", round(sum(fittedyvalues)), " Mio\n", sep=""));
          # catted data with fitted peak
          genome_size_catted <- round(sum(c(1:length(yfit2))*yfit2[1:length(yfit2)]/fitted_peak_value))
          # genomic kmers
          # cat(paste("        catted  genomic ", sizek, "mers: ", round(sum(yfit2)), " Mio\n", sep=""));
          # corrected genome size (not meaningful)
          genome_size_corrected  <- round(genome_size_catted/myxifit[1])
          # genome size: catted data + raw mean
          if(sum(hetfit[1:selected] > 0))
          {
            end_for_mean <- round(fitted_peak_value*2 - fitted_peak_value*2/4) # >7% error; k should be small
          }else
            if((sizek=="15" || sizek=="17") && myxifit[1]>=0.93 && myxifit[1]<=1.07) # normal case: reduce effect of copy-2 k-mers
            {
              end_for_mean <- round(fitted_peak_value*2 - fitted_peak_value*1/4)
            }else
              if(myxifit[1]  > 1.07)
              {
                end_for_mean <- round(fitted_peak_value*2 + fitted_peak_value*2/4)
              }else
                if(myxifit[1]  < 0.93)
                {
                  end_for_mean <- round(fitted_peak_value*2 - fitted_peak_value*4/8)
                }else
                {
                  end_for_mean <- round(fitted_peak_value*2)
                }
          #### caution: end_for_mean should be not larger than vector size
          end_for_mean <- min(end_for_mean, length(dr[,1]))       ####
          ####
          if(het_observed && only_one_het_peak) ## new ##
          {
            end_for_mean     <- round(het_peak_pos*4 - het_peak_pos*4/4)
            end_for_mean     <- min(end_for_mean, length(dr[,1])) ####
            #cat("end_for_mean (only one het peak): ", end_for_mean, '\n')
	    cat("hetfit: ceiling(abs(round(hetfit[1:end_for_mean]))/1000): ",ceiling(abs(round(hetfit[1:end_for_mean]))/1000),"\n")
            xtmp             <- rep(1:end_for_mean,ceiling(abs(round(hetfit[1:end_for_mean]))/1000))
            first_mean_raw   <- 2*mean(xtmp[xtmp>=1 & xtmp<=end_for_mean])
          }else
            if(het_observed & main_peak_is_hom==F)
            {
              het_end_for_mean <- first_peak_pos
              het_end_for_mean <- min(het_end_for_mean, length(dr[,1])) ####
              xtmp             <- rep(1:het_end_for_mean,ceiling(abs(round(hetfit[1:het_end_for_mean]))/1000))
              first_mean_raw   <- 2*mean(xtmp[xtmp>=1 && xtmp<=het_end_for_mean])
            }else
              if(het_observed && only_one_hom_peak) # caution: expert suppl.
              {
                dtmp           <- ceiling(abs(round(yfit2[1:end_for_mean]-(1-het_fitting_delta)*hetfit[1:end_for_mean]))/1000) # 2016-08-25
                ## start of specific in v1.94: from [0.5*het_peak, 1.5*het_peak]
                offset_count   <- round((fittedyvalues[1:end_for_mean]-yfit2[1:end_for_mean])/1000)
                offl <- round(0.5*which.max(hetfit))
                offr <- round(1.5*which.max(hetfit))
                for (off_i in c(offl:offr))
                {
                  if(offset_count[off_i] > 0)
                  {
                    dtmp[off_i] <- dtmp[off_i] + offset_count[off_i]
                  }
                  else
                  {
                    dtmp[off_i] <- dtmp[off_i] - offset_count[off_i]
                  }
                }
                ## end   of specific in v1.94
                xtmp           <- rep(1:end_for_mean,dtmp)
                first_mean_raw <- mean(xtmp[xtmp>=1 & xtmp<=end_for_mean])
              }else
              {
                xtmp           <- rep(1:end_for_mean,ceiling(abs(round(yfit2[1:end_for_mean]-hetfit[1:end_for_mean]))/1000))
                first_mean_raw <- mean(xtmp[xtmp>=1 & xtmp<=end_for_mean]);
              }
          #
          genome_size_corrected2 <- round(sum(c(1:length(yfit2))*yfit2[1:length(yfit2)]/first_mean_raw))
          #
          #
          ## heterozygosity est: new feature - 2016-11-18
          if(het_observed)
          {
            # alpha = 1 / [ (O/E)*K + K ]
            # het
            het_end <- round(first_mean_raw)
            hom_end <- min(round(2*first_mean_raw), selected)
            if(only_one_hom_peak)
            {
              # only for humans: do the correction for those k-mers from X-Y
              het_correction <- 0
              if(species=="human")
              {
                het_correction <- (155270560 - 59373566)*first_mean_raw/2;
              }
              #sum(c(1:het_end)*(1-het_fitting_delta)*hetfit[1:het_end])
              alpha     <- 1 / ( sum(seq(1:hom_end)*homfit[seq(100, length(homfit), 100)][1:hom_end])/  sum(seq(1:het_end)*hetfit[1:het_end])/(1-het_fitting_delta)*as.numeric(sizek) + as.numeric(sizek) )
              alpha_cor <- 0
              alpha_cor <- 1 / ( sum(seq(1:hom_end)*homfit[seq(100, length(homfit), 100)][1:hom_end])/ (sum(seq(1:het_end)*hetfit[1:het_end])-het_correction)/(1-het_fitting_delta)*as.numeric(sizek) + as.numeric(sizek) )
            }
            else
            {
              # only for humans: do the correction for those k-mers from X-Y
              het_correction <- 0
              if(species=="human")
              {
                het_correction <- (155270560 - 59373566)*first_mean_raw/2;
              }
              #sum(c(1:het_end)*hetfit[1:het_end])
              alpha     <- 1 / ( sum(seq(1:hom_end)*homfit[seq(100, length(homfit), 100)][1:hom_end])/  sum(seq(1:het_end)*hetfit[1:het_end])*as.numeric(sizek) + as.numeric(sizek) )
              alpha_cor <- 0
              alpha_cor <- 1 / ( sum(seq(1:hom_end)*homfit[seq(100, length(homfit), 100)][1:hom_end])/ (sum(seq(1:het_end)*hetfit[1:het_end])-het_correction )*as.numeric(sizek) + as.numeric(sizek) )
            }
          }
          ## end of heterozygosity rate est: new feature - 2016-11-18
          #
          #
          # single copy region of genome
          singlestart = 1
          singleend   = min(round(2*first_mean_raw),length(dr[,1])-1) #### caution: out of bounds
          singlesize  = sum(dhet[singlestart:singleend, 1]*yfit2[singlestart:singleend]/first_mean_raw)
          # begin: new on 2016-08-08
          repsize  = round(sum(dhet[(singleend+1):length(yfit2), 1]*yfit2[(singleend+1):length(yfit2)]/first_mean_raw))
          #legend(45, 1.36*max(dr[border_pos:60, 2]),
          if(het_observed)
          {
            legend("topright",
                   fill   = c('gray','white', 'white', 'cyan', 'white', 'blue',
                              'white', 'white', 'white', 'red', 'white', 'white'),
                   legend = c(paste('Observed (obs.): ',
                                    round(genome_size/10^6, digits = 3), ' Mb (error-excluded: ',
                                    round(genome_size_filtering_error/10^6, digits = 3), ' Mb)', sep=""),
                              paste('* k-mer_cov obs.: ',     first_peak_pos, sep=""),
                              paste('* signal_error_border ', het_min_valid_pos,                 sep=""),
                              paste('Heterozygous fitting (intermediate info)'),
                              paste('* het rate: ',
                                    round(alpha, digits=8),  " (ori) : (cor) ",
                                    round(alpha_cor, digits = 8),                                sep=""),
                              paste('Fitted count with fitted k-mer_cov: ',
                                    round(genome_size_fitted/10^6, digits = 3), ' Mb',           sep=""),
                              paste('* k-mer_cov fit.: ',
                                    round(fitted_peak_value, digits=2),                          sep=""),
                              paste('* 1st-sd:',
                                    round(first_sd, digits=2),                                   sep=""),
                              paste('* 1st-skewness: ',
                                    round(myxifit[1], digits=2),                                 sep=""),
                              paste('Fitted+obs. with corrected k-mer_cov: ',
                                    round(genome_size_corrected2/10^6, digits = 3), ' Mb',       sep=""),
                              paste('* k-mer_cov cor.: ',
                                    round(first_mean_raw, digits=2),                             sep=""),
                              paste('* repetitive_ratio ',
                                    round(repsize/genome_size_corrected2, digits=2), ': ',
                                    round(repsize/10^6, digits = 3), ' Mb',                      sep="")),
                   border="NA",
                   cex=0.7)
          }
          else
          {
            legend("topright",
                   fill   = c('gray','white', 'white', 'blue', 'white',
                              'white', 'white', 'red', 'white', 'white'),
                   legend = c(paste('Observed (obs.): ',
                                    round(genome_size/10^6, digits = 3), ' Mb (error-excluded: ',
                                    round(genome_size_filtering_error/10^6, digits = 3), ' Mb)',               sep=""),
                              paste('* k-mer_cov obs.: ',     first_peak_pos, sep=""),
                              paste('* signal_error_border ', border_pos,                     sep=""),
                              paste('Fitted count with fitted k-mer_cov: ',
                                    round(genome_size_fitted/10^6, digits = 3), ' Mb',                         sep=""),
                              paste('* k-mer_cov fit.: ',     round(fitted_peak_value, digits=2),              sep=""),
                              paste('* 1st-sd:',              round(first_sd, digits=2),                       sep=""),
                              paste('* 1st-skewness: ',       round(myxifit[1], digits=2),                     sep=""),
                              paste('Fitted+obs. with corrected k-mer_cov: ',
                                    round(genome_size_corrected2/10^6, digits = 3), ' Mb',                     sep=""),
                              paste('* k-mer_cov cor.: ',     round(first_mean_raw, digits=2),                 sep=""),
                              paste('* repetitive_ratio ',
                                    round(repsize/genome_size_corrected2, digits=2), ': ',
                                    round(repsize/10^6, digits = 3), ' Mb',                                    sep="")),
                   border="NA",
                   cex=0.7)
          }
          # end: new on 2016-08-08
          # genome sizes
          size_all[i]  <- genome_size
          size_exl[i]  <- genome_size_filtering_error
          size_fit[i]  <- genome_size_fitted
          size_cat[i]  <- genome_size_catted
          size_cor[i]  <- genome_size_corrected
          size_cor2[i] <- genome_size_corrected2
          # rbind(size_all, size_exl, size_fit, size_cat)
          # genome-wide average kmer coverage related info
          kmer_peak_cnt[i]  <- dr[first_peak_pos, 2]
          kmer_peak_pos[i]  <- first_peak_pos
          kmer_peak_pos2[i] <- fitted_peak_value
        }
        # labels for later plotting purpose
        labels[3*i-1]  <- paste("", sizek, sep="")
        namebank[i]    <- sizek
      }
      else
      {
        dev.off()
        stop(paste("file ", countingfile, " not found!",sep=""))
        quit("no")
      }
    }
    # summarize all kmer counting plots
    ylimmax <- 0
    for (sizek in targetsizek[1]) {
      histx<-read.table(countingfile)
      peaks2             <- findpeaks(histx$V2[3:selected])
      peaks2[,2:4]       <- peaks2[,2:4]+2
      if(1.5*max(peaks2[,1]) > ylimmax)
      {
        ylimmax <- 1.5*max(peaks2[,1])
      }
    }
    # colors
    colors <- rainbow(length(size_all))
    #
    ## plot comparision among gsize-estimations using different info
    genome_size_summary  <-rbind(size_all, size_exl, size_cat, size_fit, size_cor2)
    genome_size_summary2 <-rbind(size_exl, size_fit, size_cor2)
    cat(paste("Genome size estimate for ", sample, ": ", size_cor2, " bp.\n\n", sep=""))
    # caution
    if(length(which(is.na(genome_size_summary))) > 0)
    {
      cat("\n caution: some samples have NA predicted for observed genome size!\n")
      genome_size_summary[is.na(genome_size_summary)]   <- 0
    }
    if(length(which(is.na(genome_size_summary2))) > 0)
    {
      cat("\n caution: some samples have NA predicted for fitted   genome size!\n")
      genome_size_summary2[is.na(genome_size_summary2)] <- 0
    }
    if(length(genome_size_summary[1,])!=0 & length(genome_size_summary2[1,])!=0)
    {
      # record in file
      write.table(genome_size_summary,
                  paste(path, '/',prefix, vers, 'est.', sample, ".genome.size.estimated.k",
                        min(targetsizek), 'to', max(targetsizek), ".fitted.txt", sep=""),
                  quote = FALSE, row.names = TRUE, col.names=FALSE)
      # bar plot
      mingsize <- min(genome_size_summary2[!is.na(genome_size_summary2)])
      maxgsize <- max(genome_size_summary2[!is.na(genome_size_summary2)])
      colors <- grey((2:0)/3)
      mp <- barplot(genome_size_summary2, beside=TRUE, col=colors,
                    xlab=paste("Size k", sep=""), ylab="Genome size (bp)",
                    axes = FALSE, axisnames  = FALSE, border=NA,
                    ylim=c(0.5*max(genome_size_summary2),1.2*max(genome_size_summary2)),
                    xpd=FALSE)
      legend("topright",
             legend = c(paste('EXL: ', round(mean(genome_size_summary2[1,]/10^6), digits=3), 'Mb'),
                        paste('FIT: ', round(mean(genome_size_summary2[2,]/10^6), digits=3), 'Mb'),
                        paste('COR: ', round(mean(genome_size_summary2[3,]/10^6), digits=3), 'Mb')),
             fill = colors, cex=.8, border=NA)
      text(mp, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.8)
      axis(2, cex.axis=.8)
      title(main=paste("Genome size estimation by error-excluding, fitting, and correcting\n", sep=""))
      dev.off()
    }
  } # end of sample
  # output est. on het if asked.
  if(het_observed)
  {
    write(paste("Het_rate ", round(alpha, digits = 8), " ", round(alpha_cor, digits = 8), sep=""),
          file = paste(path, '/',prefix, vers, 'est.',
                       sample, ".genome.size.estimated.k",
                       min(targetsizek), 'to', max(targetsizek),
                       ".fitted.txt", sep=""),
          append = T,
          sep = " ")
  }
  # output est. on ratio of repeats
  write(paste("Est. ratio of repeats ", round(repsize/genome_size_corrected2, digits=8), sep=""),
        file = paste(path, '/',prefix, vers, 'est.',
                     sample, ".genome.size.estimated.k",
                     min(targetsizek), 'to', max(targetsizek),
                     ".fitted.txt", sep=""),
        append = T,
        sep = " ")
  # output est. on avg k-mer cov
  write(paste("Final k-mer cov ", round(first_mean_raw, digits=8), sep=""),
        file = paste(path, '/',prefix, vers, 'est.',
                     sample, ".genome.size.estimated.k",
                     min(targetsizek), 'to', max(targetsizek),
                     ".fitted.txt", sep=""),
        append = T,
        sep = " ")
  if(het_observed && only_one_hom_peak)
  {
    write(paste("het_fitting_delta ", round(het_fitting_delta, digits=8), sep=""),
          file = paste(path, '/',prefix, vers, 'est.',
                       sample, ".genome.size.estimated.k",
                       min(targetsizek), 'to', max(targetsizek),
                       ".fitted.txt", sep=""),
          append = T,
          sep = " ")
  }
  #
  end_t <- Sys.time()
  cat('Time consumed: ', format(end_t-start_t), '\n')
  return(list(het_xfit_right_save,hom_peak_pos_save,het_peak_pos_save))
}
##
## additional
filter_peaks <- function(peaks, histo)
{
  # find out a major peak with enough support information on both sides (to exlcude potential local maxima)
  pos <- 1 # peak position (namely, the k-mer freq at the peak)
  i   <- 0 # traverse candidate peaks
  # a peak freq=pos should start at least from 5, or else too small to perform curve fitting
  while(pos < 5 & i < length(peaks[,1]))
  {
    i <- i+1
    pos <- peaks[i, 2] # x-axis peak
  }
  if(pos < 5)
  {
    stop(paste("   Error: k-mer coverage seemed too low to perform fitting. Minimum average k-mer coverage 5 required. Program exited.\n",sep=""))
  }
  # check if there are enough support on both sides of pos
  while(i < length(peaks[,1]))
  {
    # check left: how many on left are smaller than count at current peak
    posleft <- pos - 1
    j <- 0 # counter of smaller cases on left
    while(posleft>=1)
    {
      if(histo[posleft, 2] < histo[pos, 2])
      {
        j <- j + 1;
      }
      posleft <- posleft - 1
    }
    # check right: if there are right-side counts larger than count at the current peak
    posright <- pos + 1
    k <- 0 # counter of larger cases on right
    maxpos <- peaks[length(peaks[,1]), 2]
    while(posright < maxpos)
    {
      if(histo[posright, 2] > histo[pos, 2])
      {
        k <- k + 1
        if(k>4)
        {
          break
        }
      }
      posright <- posright + 1
    }
    # condition to continue or break
    if(j <= 4 || k>4)
    {
      i <- i + 1
      pos <- peaks[i, 2] # x-axis peak
    }
    else
    {
      break
    }
  }
  # filter out peaks near erroneous peak
  rpeaks <- peaks
  fpi <- 1
  while(rpeaks[fpi, 2] < pos/4)
  {
    rpeaks <- rpeaks[-fpi,]
  }
  return(rpeaks)
}

get_het_pos <- function(histo_data){
        data <- histo_data
        y_data <- data$V2

        # Calculate slope (second derivative)
        slope <- diff(y_data)

        # Find the index where slope becomes positive and its following 5 slopes remain positive or zero
        start_index <- NULL
        for (i in 2:(length(slope) - 4)) {
          if (slope[i] > 0 && all(slope[(i+1):(i+4)] >= 0)) {
            start_index <- i
            break
          }
        }
        cat("start index: ",start_index,"\n")
        end_index <- NULL
        # Find the first interval where slope starts decreasing after positive slope
        for (i in (start_index + 1):(length(slope) - 1)) {
          if (slope[i] < 0 && slope[i+1] < 0) {
            end_index <- i
            break
          }
          else if (slope[i] >= 0 && slope[i] < slope[i-1] && slope[i] < slope[i+1] &&
                        slope[i-1] < slope[i-2] && slope[i-2] < slope[i-3] &&
              slope[i+1] < slope[i+2] && slope[i+2] < slope[i+3]){
            if (is.null(end_index)) {
              end_index <- i
              break
            }
          }
        }

        #print(end_index)
        return(list(start_index ,end_index))
}
# => *_21mer.histo.genome.size.estimated.k21to21.fitted_hetfit_count.txt needed by GSE_sp.R

##  merge function
#path="/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a2_initial_assembly/"
#path="/home/biodata/group_sun/reads/proj_solanum_tuberosum/wgs_4_longshu/"

################################################## main #############################################################
#' @title Estimate genome size of polyploid species using k-mer frequencies.
#'
#' @description findGSEX is a function for (up to 10 ploidy)
#' genome size estimation by fitting k-mer frequencies iteratively
#' with a normal distribution model. (version still under testing)
#'
#' @description To use findGSEX, one needs to prepare a histo file,
#' which contains two tab-separated columns.
#' The first column gives frequencies at which k-mers occur in reads,
#' while the second column gives counts of such distinct k-mers.
#' Parameters k and related histo file are required for any estimation.
#'
#' @description Dependencies (R library) required: pracma, fGarch - see INSTALL.
#'
#'
#' @param histo is the histo file (mandatory).
#' @param sizek is the size of k used to generate the histo file (mandatory).
#' K is involved in calculating heterzygosity if the genome is heterozygous.
#' @param outdir is the path to write output files (optional).
#' If not provided, by default results will be written in the folder
#' where the histo file is.
#' @param exp_hom a rough average k-mer coverage for finding the homozygous regions.
#' In general, one can get peaks in the k-mer frequencies file, but has to
#' determine which one is for the homozygous regions, and which one is for the
#' heterozygous regions. It is optional, however, it must be provided
#' if one wants to estimate size for a heterozygous genome.
#' VALUE for exp_hom must satisfy fp < VALUE < 2*fp, where fp is the freq for homozygous peak.
#' If not provided, 0 by default assumes the genome is homozygous.
#' @export
#'
findGSEX <- function(path, samples, sizek, exp_hom, ploidy, range_left, range_right, xlimit, ylimit ,output_dir="outfile"){

  if (!grepl("/$", path)) {
    path <- paste0(path, "/")
  }
  if (!grepl("/$", output_dir)) {
    output_dir <- paste0(output_dir, "/")
  }

  if(missing(exp_hom)) exp_hom <- 0
  if(missing(ploidy)) ploidy <- 2

  if(ploidy > 7) {ploidy_assign <- FALSE
  }
  else{
    ploidy_assign <- TRUE}
  if(missing(range_left)) range_left <- 6
  if(missing(range_right)) range_right <- 9
  if(missing(sizek)) sizek <- 21

  #cat(range_left,'\n')
  if(missing(output_dir))   output_dir  <- getwd()

  if(missing(samples))
  {
    stop(paste("Please provide histo file name, ", " Program exited.\n",sep=""))
  }
  if(missing(sizek))
  {
    stop(paste("Cannot find sizek info: ", sizek, ". Program exited.\n",sep=""))
  }
  ploidy_raw <- 0
  for (sample in samples) {

    pdf(paste(path, output_dir ,sample,"_hap_genome_size_est.pdf", sep=""), family="Helvetica", height=4, width=7.08661) # need to update

    if(ploidy <= 2){
      cat('Ploidy less than 2 ,starting...\n')
      brewer_palette <- brewer.pal(ploidy, "Set1")
      findGSE_raw(histo=paste0(path ,sample),sizek=sizek, outdir=paste0(path, output_dir), exp_hom=exp_hom)

      cat('Start findGSE plot drawing, please wait...\n')
      ## read data

      histo_raw <- read.table(paste0(path ,sample))
      cat('histo_raw read done...\n')
      histo_fit <- read.table(paste(path, output_dir,'v1.95.est.',
                          sample, ".genome.size.estimated.k",
                          min(sizek), 'to', max(sizek),".fitted_fullfit_count.txt",sep=""))

      cat('histo_fit read done...\n')
      if(exp_hom != 0){
        histo_het <- read.table(paste(path, output_dir,'v1.95.est.',
                          sample, ".genome.size.estimated.k",
                          min(sizek), 'to', max(sizek),".fitted_hetfit_count.txt",sep=""))
        cat('histo_het read done...\n')
      }
      else{
        histo_het <- read.table(paste(path, output_dir,'v1.95.est.',
                          sample, ".genome.size.estimated.k",
                          min(sizek), 'to', max(sizek),".fitted_hetfit_count.txt",sep=""))
        histo_first_fit_hom <-   histo_het
        # histo_first_fit_hom <- read.table(paste(path, output_dir,'v1.95.est.',
        #                   sample, ".genome.size.estimated.k",
        #                   min(sizek), 'to', max(sizek),".fitted_first_fit_hom_count.txt",sep=""))
        cat('first_fit_hom read done...\n')
      }


      #load data
      load("variable_data.Rdata")

      #
      genome_size <- data_to_save$genome_size
      genome_size_filtering_error <- data_to_save$genome_size_filtering_error
      first_peak_pos <- data_to_save$first_peak_pos
      genome_size_fitted <- data_to_save$genome_size_fitted
      fitted_peak_value <- data_to_save$fitted_peak_value
      first_sd <- data_to_save$first_sd
      myxifit <- data_to_save$myxifit
      genome_size_corrected2 <- data_to_save$genome_size_corrected2
      first_mean_raw <- data_to_save$first_mean_raw
      repsize <- data_to_save$repsize

      ## plot
      #plot
      if(xlimit == -1 ){
        xlimit = first_peak_pos * ploidy * 2
      }
      if(ylimit == -1){
        ylimit = max(histo_fit) * 1.1
      }
      plot(histo_raw$V1,
          histo_raw$V2,
          xlim=c(1, xlimit),
          ylim=c(1, ylimit),
          xlab="Coverage of k-mers",
          ylab="Frequency",
          type="b",
          col ="gray",
          lwd=1,
          cex=1.5,
          frame.plot = F)
      polygon(c(histo_raw$V1, rev(histo_raw$V1)),
      c(histo_raw$V2, rep(1, length(histo_raw$V2))),
      col = alpha("gray", 0.5), border = NA)
      # final used for gse
      lines(histo_fit$V1,  histo_fit$V2, col=alpha("orangered", 0.8), lwd=6)
      if(exp_hom != 0){
      # het fitting
      lines(histo_het$V1,  histo_het$V2, col="cyan", lwd=3, lty=2)

      polygon(c(histo_het$V1, rev(histo_het$V1)),
      c(histo_het$V2, rep(1, length(histo_het$V2))),
      col = alpha("cyan", 0.3), border = NA)
      }
      cat('plot het line done...\n')

      polygon(c(histo_fit$V1, rev(histo_fit$V1)),
      c(histo_fit$V2, rep(1, length(histo_fit$V2))),
      col = alpha("orangered", 0.3), border = NA)






      #
      if (exp_hom != 0) {
        het_min_valid_pos <- data_to_save$het_min_valid_pos
        alpha <- data_to_save$alpha
        alpha_cor <- data_to_save$alpha_cor

        #
      } else {
        border_pos <- data_to_save$border_pos
      }

      ## add line 8031 add hom peak
      hom_uplimit <- round(first_mean_raw * 2)

      if(exp_hom != 0){
          hom_histo_unique_region <- histo_fit[1:hom_uplimit,2] - histo_het[1:hom_uplimit,2]
          hom_histo_unique_region <- as.data.frame(cbind(1:hom_uplimit,hom_histo_unique_region))
          colnames(hom_histo_unique_region) <- c("V1","V2")
          first_mean_raw_half <- round(first_mean_raw / 2)
          lines(hom_histo_unique_region[first_mean_raw_half:hom_uplimit,1],  hom_histo_unique_region[first_mean_raw_half:hom_uplimit,2], col=brewer_palette[2], lwd=3, lty=2)
      }
      else{
          hom_histo_unique_region <- histo_first_fit_hom
          # hom_histo_unique_region <- as.data.frame(cbind(1:hom_uplimit,hom_histo_unique_region))
          colnames(hom_histo_unique_region) <- c("V1","V2")
          first_mean_raw_half <- round(first_mean_raw / 2)
          lines(hom_histo_unique_region[1:round(first_mean_raw * 2),1],  hom_histo_unique_region[1:round(first_mean_raw * 2),2], col=brewer_palette[2], lwd=3, lty=2)

      }
      # first_mean_raw_half <- round(first_mean_raw / 2)
      # lines(hom_histo_unique_region[first_mean_raw_half:hom_uplimit,1],  hom_histo_unique_region[first_mean_raw_half:hom_uplimit,2], col=brewer_palette[2], lwd=3, lty=2)
      genome_size_hom <- sum(hom_histo_unique_region$V1 * hom_histo_unique_region$V2/1e+06)/ first_mean_raw
      genome_size_hom <- round(genome_size_hom, digits = 3)
      if(exp_hom != 0)
          {
            legend("topright",
                   fill   = c('gray','white', 'white', 'cyan','white',brewer_palette[2],'white',
                               "orangered",'white','white'),
                   legend = c(paste('Observed (obs.): ',
                                    round(genome_size/10^6, digits = 3), ' Mb (error-excluded: ',
                                    round(genome_size_filtering_error/10^6, digits = 3), ' Mb)', sep=""),
                              paste('* k-mer_cov obs.: ',     first_peak_pos, sep=""),
                              paste('* signal_error_border ', het_min_valid_pos,                 sep=""),
                              paste('Heterozygous fitting (intermediate info)'),
                              paste('* het rate: ',
                                    round(alpha, digits=8),  " (ori) : (cor) ",
                                    round(alpha_cor, digits = 8),                                sep=""),
                              paste('Unique hom-region: ',
                                    genome_size_hom, ' Mb',       sep=""),
                              paste('Fitted+obs. with corrected k-mer_cov: ',
                                    round(genome_size_corrected2/10^6, digits = 3), ' Mb',       sep=""),
                              paste('* k-mer_cov cor.: ',
                                    round(first_mean_raw, digits=2),                             sep=""),
                              paste('* repetitive_ratio ',
                                    round(repsize/genome_size_corrected2, digits=2), ': ',
                                    round(repsize/10^6, digits = 3), ' Mb',                      sep="")),
                   border="NA",
                   cex=0.5)
          }
          else
          {
            legend("topright",
                   fill   = c('gray','white', 'white',brewer_palette[2],'white', "orangered",'white','white'),
                   legend = c(paste('Observed (obs.): ',
                                    round(genome_size/10^6, digits = 3), ' Mb (error-excluded: ',
                                    round(genome_size_filtering_error/10^6, digits = 3), ' Mb)',               sep=""),
                              paste('* k-mer_cov obs.: ',     first_peak_pos, sep=""),
                              paste('* signal_error_border ', border_pos,                     sep=""),
                              paste('Unique hom-region: ',
                                    genome_size_hom, ' Mb',       sep=""),
                              paste('Fitted+obs. with corrected k-mer_cov: ',
                                    round(genome_size_corrected2/10^6, digits = 3), ' Mb',                     sep=""),
                              paste('* k-mer_cov cor.: ',     round(first_mean_raw, digits=2),                 sep=""),
                              paste('* repetitive_ratio ',
                                    round(repsize/genome_size_corrected2, digits=2), ': ',
                                    round(repsize/10^6, digits = 3), ' Mb',                                    sep="")),
                   border="NA",
                   cex=0.5)
          }

      dev.off()

      end_time <- proc.time()
      execution_time <- end_time - start_time
      print(paste("Program running time:", sum(execution_time)))
      ## draw png file
      png(paste(path, output_dir ,sample,"_hap_genome_size_est.png", sep=""),  width = 1200, height = 800, res = 200)
      par(family = "Helvetica")

      plot(histo_raw$V1,
          histo_raw$V2,
          xlim=c(1, xlimit),
          ylim=c(1, ylimit),
          xlab="Coverage of k-mers",
          ylab="Frequency",
          type="b",
          col ="gray",
          lwd=1,
          cex=1.5,
          frame.plot = FALSE,
          cex.axis = 0.6)
      polygon(c(histo_raw$V1, rev(histo_raw$V1)),
      c(histo_raw$V2, rep(1, length(histo_raw$V2))),
      col = alpha("gray", 0.5), border = NA)
      # final used for gse
      lines(histo_fit$V1,  histo_fit$V2, col=alpha("orangered", 0.8), lwd=6)
      if(exp_hom != 0){

          lines(hom_histo_unique_region[first_mean_raw_half:hom_uplimit,1],  hom_histo_unique_region[first_mean_raw_half:hom_uplimit,2], col=brewer_palette[2], lwd=3, lty=2)
      }
      else{

          lines(hom_histo_unique_region[1:round(first_mean_raw * 2),1],  hom_histo_unique_region[1:round(first_mean_raw * 2),2], col=brewer_palette[2], lwd=3, lty=2)

      }

      if(exp_hom != 0){
      # het fitting
      lines(histo_het$V1,  histo_het$V2, col="cyan", lwd=3, lty=2)

      polygon(c(histo_het$V1, rev(histo_het$V1)),
      c(histo_het$V2, rep(1, length(histo_het$V2))),
      col = alpha("cyan", 0.3), border = NA)
      }

      polygon(c(histo_fit$V1, rev(histo_fit$V1)),
      c(histo_fit$V2, rep(1, length(histo_fit$V2))),
      col = alpha("orangered", 0.3), border = NA)


      if(exp_hom != 0)
          {
            legend("topright",
                   fill   = c('gray','white', 'white', 'cyan','white',brewer_palette[2],'white',
                               "orangered",'white','white'),
                   legend = c(paste('Observed (obs.): ',
                                    round(genome_size/10^6, digits = 3), ' Mb (error-excluded: ',
                                    round(genome_size_filtering_error/10^6, digits = 3), ' Mb)', sep=""),
                              paste('* k-mer_cov obs.: ',     first_peak_pos, sep=""),
                              paste('* signal_error_border ', het_min_valid_pos,                 sep=""),
                              paste('Heterozygous fitting (intermediate info)'),
                              paste('* het rate: ',
                                    round(alpha, digits=8),  " (ori) : (cor) ",
                                    round(alpha_cor, digits = 8),                                sep=""),
                              paste('Unique hom-region: ',
                                    genome_size_hom, ' Mb',       sep=""),
                              paste('Fitted+obs. with corrected k-mer_cov: ',
                                    round(genome_size_corrected2/10^6, digits = 3), ' Mb',       sep=""),
                              paste('* k-mer_cov cor.: ',
                                    round(first_mean_raw, digits=2),                             sep=""),
                              paste('* repetitive_ratio ',
                                    round(repsize/genome_size_corrected2, digits=2), ': ',
                                    round(repsize/10^6, digits = 3), ' Mb',                      sep="")),
                   border="NA",
                   cex=0.5)
          }
          else
          {
            legend("topright",
                   fill   = c('gray','white', 'white',brewer_palette[2],'white', "orangered",'white','white'),
                   legend = c(paste('Observed (obs.): ',
                                    round(genome_size/10^6, digits = 3), ' Mb (error-excluded: ',
                                    round(genome_size_filtering_error/10^6, digits = 3), ' Mb)',               sep=""),
                              paste('* k-mer_cov obs.: ',     first_peak_pos, sep=""),
                              paste('* signal_error_border ', border_pos,                     sep=""),
                              paste('Unique hom-region: ',
                                    genome_size_hom, ' Mb',       sep=""),
                              paste('Fitted+obs. with corrected k-mer_cov: ',
                                    round(genome_size_corrected2/10^6, digits = 3), ' Mb',                     sep=""),
                              paste('* k-mer_cov cor.: ',     round(first_mean_raw, digits=2),                 sep=""),
                              paste('* repetitive_ratio ',
                                    round(repsize/genome_size_corrected2, digits=2), ': ',
                                    round(repsize/10^6, digits = 3), ' Mb',                                    sep="")),
                   border="NA",
                   cex=0.5)
          }

      dev.off()

    }
    else{
        cat("ploidy id larger that 2, enter next process...",path,sample,"\n")
	      histo_data <- read.table(paste(path, sample, sep = ""), header = FALSE)
        het_pos_list <- get_het_pos(histo_data)
        start_pos <- het_pos_list[[1]]
        het_pos <- het_pos_list[[2]]
        # left_fit_ratio_list <- c(0.835, 0.855, 0.825)
        left_fit_ratio_list <- c(0.835)
        ## define selected color
        # select color
        brewer_palette <- brewer.pal(ploidy, "Set1")
        if(het_pos < 80){
            cat("het_pos is ",het_pos," enter het_pos < 80 process... \n")
            portion_size <- c()
            success_file <- TRUE
            continue_reverse_flag = FALSE
            for (ploidy_ind in seq(ploidy)){
            if(ploidy_ind == 1){
              exp_hom_temp = exp_hom
              avg_cov = 0
              cat('ploidy index : ',ploidy_ind)
              histo_raw <- read.table(paste(path, sample, sep="") ) # from jellyfish
              ## 1204 add cutting
              if((het_pos - start_pos) < 8){
                histo_org <- histo_raw
                histo_org[1:(start_pos-1),2] <- 1
                # pivot_value <- histo_org$V2[start_pos]
                # histo_org$V2[1:(start_pos-1)] <- pivot_value - (histo_org$V2[1:(start_pos-1)] - pivot_value)


                # histo_org[(het_pos-29+het_pos):het_pos, 2] <- rev(histo_org[round(het_pos):29, 2])

                # histo_org[30:2000, 2] <- 1
                write.table(histo_org, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = F)
              }
              ##end 1204
              for(left_fit_ratio in left_fit_ratio_list){
                   if((het_pos - start_pos) < 8){
                      cat("read cutting file...\n")
                      fit_para_save_list = findGSE_sp(histo = paste(path, sample,"_",ploidy_ind, sep = ""), sizek = sizek , outdir = paste(path, output_dir, sep = ""), exp_hom = exp_hom_temp,ploidy_ind=ploidy_ind, avg_cov=avg_cov, left_fit_ratio=left_fit_ratio)
                      het_xfit_right_save <- fit_para_save_list[[1]]
                      hom_peak_pos_save <- fit_para_save_list[[2]]
                      het_peak_pos_save <- fit_para_save_list[[3]]
                      histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample, "_", ploidy_ind, ".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R
                   }
                  else{
                    fit_para_save_list = findGSE_sp(histo = paste(path, sample, sep = ""), sizek = sizek , outdir = paste(path, output_dir, sep = ""), exp_hom = exp_hom_temp,ploidy_ind=ploidy_ind, avg_cov=avg_cov, left_fit_ratio=left_fit_ratio)
                    het_xfit_right_save <- fit_para_save_list[[1]]
                    hom_peak_pos_save <- fit_para_save_list[[2]]
                    het_peak_pos_save <- fit_para_save_list[[3]]
                    histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R

                  }                  #png(paste("apple_hap_genome_size_est.jpg", sep=""),  height=400, width=708)
                  ####
                                    #print(histo_raw[1:100, ])#
                  #print(histo_het[1:100, ])
                  len       <- length(histo_raw$V1)
                  # find the haplotype average kmer coverage.
                  happeak_index <- which.max(histo_het$V2) # coverage at het-peak 14, x, ...
                  #cat("hap peak at", happeak_index, '\n')
                  #cat("left peak at", range_left, '\n')
                  range_avg_cov <- (happeak_index-range_left):(happeak_index+range_right)
                  avg_cov       <- round(sum(histo_raw[range_avg_cov, 1] * histo_raw[range_avg_cov, 2] / sum(histo_raw[range_avg_cov, 2])), digits = 2)
                  cat('avg_cov: ',avg_cov,'\n')
                  final_fit_peak_pos_ind <- which.max(histo_het$V2)
                  final_fit_peak_pos <-  histo_het[final_fit_peak_pos_ind,1]
                  distance_fit_two_line <- abs(avg_cov*ploidy_ind - final_fit_peak_pos)
                  cat('difference of cyan peaks: ',final_fit_peak_pos,', gray peak: ',avg_cov*ploidy_ind,'\n')
                  if(distance_fit_two_line <= 5){
                      cat(left_fit_ratio,' meet the requirement, continue...\n')
                      break
                  }
              }



              ## estimate the size of haploid portion
              hap_size <- sum(as.numeric(histo_het$V1) * as.numeric(histo_het$V2) / 1000000) / (avg_cov * ploidy)
              cat('portion_size for ploidy_ind 1: ',hap_size,'\n')
              portion_size <- hap_size

              # find the valley on the left of het region
              valley_index  <- which.min(histo_raw[1:happeak_index, 2])
              #cat(valley_index, '\n')
              #
              # make new histogram with fittings at het-kmers
              histo_fit                    <- histo_raw
              histo_fit[1:valley_index, 2] <- histo_het[1:valley_index, 2]
              #
              if((het_pos - start_pos) < 7){
                cat('het_pos - start_pos) < 7 ohhhhh myyyy!!\n')
                histo_fit[1:het_pos, 2] <- histo_het[1:het_pos, 2]
              }

              #plot
              if(xlimit == -1 ){
                 xlimit = avg_cov * ploidy * 2
              }
              if(ylimit == -1){
                ylimit = max(histo_fit) * 1.1
              }
              plot(histo_raw$V1,
                  histo_raw$V2,
                  xlim=c(1, xlimit),
                  ylim=c(1, ylimit),
                  xlab="Coverage of k-mers",
                  ylab="Frequency",
                  type="b",
                  col ="gray",
                  lwd=1,
                  cex=1.5,
                  frame.plot = F)
              polygon(c(histo_raw$V1, rev(histo_raw$V1)),
              c(histo_raw$V2, rep(1, length(histo_raw$V2))),
              col = alpha("gray", 0.5), border = NA)
              # final used for gse
              lines(histo_fit$V1,  histo_fit$V2, col=alpha("orangered", 0.8), lwd=6)
              # het fitting
              lines(histo_het$V1,  histo_het$V2, col="cyan", lwd=3, lty=2)

              polygon(c(histo_fit$V1, rev(histo_fit$V1)),
              c(histo_fit$V2, rep(1, length(histo_fit$V2))),
              col = alpha("orangered", 0.3), border = NA)


              polygon(c(histo_het$V1, rev(histo_het$V1)),
              c(histo_het$V2, rep(1, length(histo_het$V2))),
              col = alpha("cyan", 0.3), border = NA)

              final_fit_peak_pos_ind <- which.max(histo_het$V2)
              final_fit_peak_pos <-  histo_het[final_fit_peak_pos_ind,1]
              cat('difference of cyan peaks: ',final_fit_peak_pos,', gray peak: ',avg_cov,'\n')

              #
              abline(v = avg_cov*1, lty=3, col="gray")

              #
              genome_size_raw = sum(histo_raw[valley_index:len, 1] * histo_raw[valley_index:len, 2] / 1000000) / (ploidy*avg_cov)
              genome_size     = sum(histo_fit$V1 * (histo_fit$V2 / 1000000)) / (ploidy*avg_cov) # tetraploid size: 3374.929 Mb, haploid size: 843.7321
              #

              #
              if(ploidy == 1){
                dev.off()
                break
              }



            }
            else{
              cat("ploidy_ind is: ",ploidy_ind)
              if(!success_file){
                next
              }
              ## update histo_file
              histo_tmp = histo_raw
              het_length = dim(histo_het)[1]
              histo_tmp[1:het_length,2] = abs(histo_tmp[1:het_length,2] - histo_het[,2])
        #histo_tmp[1:(ploidy_ind/(ploidy_ind+1)*avg_cov*ploidy_ind)*0.6,2] = 0

              ## add reverse data part


              target_start_index <- round(hom_peak_pos_save)
              target_end_index <- round(hom_peak_pos_save) * 2 - round(het_peak_pos_save)
              if(ploidy_ind > 2){
                  max_het_ind <- which(histo_tmp[,2] == max(histo_tmp[round(hom_peak_pos_save):length(histo_tmp$V1),2]))
                  cat('max_het_ind is : ',max_het_ind,'\n')
                  if (max_het_ind && round(hom_peak_pos_save) < max_het_ind && (max_het_ind - round(hom_peak_pos_save)) > 3 ){
                      cat('match latter peak pattern, reverse data in ... ploidy ',ploidy_ind,'\n')
                      histo_tmp_raw <- histo_tmp
                      histo_tmp_raw[target_start_index:length(histo_tmp$V1),2] <- 1
                      histo_tmp_raw[target_start_index:target_end_index, 2] <- rev(histo_tmp_raw[round(het_peak_pos_save):round(hom_peak_pos_save), 2])
                      write.table(histo_tmp_raw, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = F)
                      #lines(histo_tmp_raw$V1,  histo_tmp_raw$V2, col="black", lwd=3)
                      continue_reverse_flag = TRUE
                  }
                  else{
                    if(continue_reverse_flag){
                    histo_tmp_raw <- histo_tmp
                    # histo_tmp_raw[which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))):length(histo_tmp_raw$V1),2] <- 1
                    histo_tmp_raw[target_start_index:length(histo_tmp$V1),2] <- 1
                    histo_tmp_raw[target_start_index:target_end_index, 2] <- rev(histo_tmp_raw[round(het_peak_pos_save):round(hom_peak_pos_save), 2])

                    #cat('pppploidy_ind range: ',which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))),'\n')
                    write.table(histo_tmp_raw, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = F)

                    }
                    else{
                    histo_tmp_raw <- histo_tmp
                    cat('This difference is ',histo_tmp[avg_cov*ploidy_ind,2] / histo_tmp[het_xfit_right_save-1,2],'\n')
                    cat('The two peak difference is ',histo_tmp[avg_cov*(ploidy_ind+1),2] / histo_tmp[avg_cov*ploidy_ind,2],'\n')
                    cat('This peak difference is ',histo_tmp[avg_cov*(ploidy_ind+1),2] / histo_tmp[het_xfit_right_save-1,2],'\n')
                    if(histo_tmp_raw[avg_cov*ploidy_ind,2] / histo_tmp_raw[het_xfit_right_save-1,2] < 3 && histo_tmp_raw[avg_cov*ploidy_ind,2] / histo_tmp_raw[het_xfit_right_save-1,2] > 1.6){
                        cat('Difference too close, lower the histo_tmp file for ',ploidy_ind,'\n')
                        histo_tmp_raw[1:het_xfit_right_save,2] <- 1
                    #histo_tmp[avg_cov*ploidy_ind:het_xfit_right_save,2] <- round(histo_tmp[avg_cov*ploidy_ind,2])
                    }
                    # histo_tmp_raw[which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))):length(histo_tmp_raw$V1),2] <- 1
                    #histo_tmp_raw[target_start_index:length(histo_tmp$V1),2] <- 1
                    #histo_tmp_raw[target_start_index:target_end_index, 2] <- rev(histo_tmp_raw[round(het_peak_pos_save):round(hom_peak_pos_save), 2])

                    #cat('pppploidy_ind range: ',which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))),'\n')
                    write.table(histo_tmp_raw, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = F)
                    # lines(histo_tmp_raw$V1,  histo_tmp_raw$V2, col="black", lwd=3)
                    }

                  }

              }
              else{
                  cat('This difference is ',histo_tmp[avg_cov*ploidy_ind,2] / histo_tmp[het_xfit_right_save-1,2],'\n')
                  cat('The two peak difference is ',histo_tmp[avg_cov*(ploidy_ind+1),2] / histo_tmp[avg_cov*ploidy_ind,2],'\n')
                  cat('This peak difference is ',histo_tmp[avg_cov*(ploidy_ind+1),2] / histo_tmp[het_xfit_right_save-1,2],'\n')
                  if(histo_tmp[avg_cov*ploidy_ind,2] / histo_tmp[het_xfit_right_save-1,2] < 3 && histo_tmp[avg_cov*ploidy_ind,2] / histo_tmp[het_xfit_right_save-1,2] > 1.6){
                    cat('Difference too close, lower the histo_tmp file for ',ploidy_ind,'\n')
                    histo_tmp[1:het_xfit_right_save,2] <- 1
                    #histo_tmp[avg_cov*ploidy_ind:het_xfit_right_save,2] <- round(histo_tmp[avg_cov*ploidy_ind,2])
                  }
                  # lines(histo_tmp$V1,  histo_tmp$V2, col="black", lwd=3)
                  write.table(histo_tmp, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = F)
              }
              ## end add reverse data part





              ## may delete
              #write.table(histo_tmp, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = F)

              exp_hom_temp = avg_cov * (ploidy_ind+1) + avg_cov * 0.1
              for(left_fit_ratio in left_fit_ratio_list){


                  fit_para_save_list = findGSE_sp(histo = paste(path, sample,"_",ploidy_ind, sep = ""), sizek = sizek , outdir = paste(path, output_dir, sep = ""), exp_hom = exp_hom_temp,ploidy_ind=ploidy_ind,avg_cov=avg_cov,left_fit_ratio=left_fit_ratio)
                  het_xfit_right_save <- fit_para_save_list[[1]]
                  hom_peak_pos_save <- fit_para_save_list[[2]]
                  het_peak_pos_save <- fit_para_save_list[[3]]
                  histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,"_",ploidy_ind,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R

                  all_slope_data <- abs(diff(diff(histo_het$V2))[(round(avg_cov)*ploidy_ind-10):(round(avg_cov)*ploidy_ind+10)])
                  if(!all(all_slope_data > 100)){
                    cat('index is :',ploidy_ind,' slope is :',all_slope_data,'\n')
                    success_file <- FALSE
                    stoped_ploidy_ind = ploidy_ind - 1
                    break
                  }

                  final_fit_peak_pos_ind <- which.max(histo_het$V2)
                  final_fit_peak_pos <-  histo_het[final_fit_peak_pos_ind,1]
                  distance_fit_two_line <- abs(avg_cov*ploidy_ind - final_fit_peak_pos)
                  cat('difference of cyan peaks: ',final_fit_peak_pos,', gray peak: ',avg_cov*ploidy_ind,'\n')
                  if(distance_fit_two_line <= 5){
                      cat(left_fit_ratio,' meet the requirement, continue...\n')
                      break
                  }

              }
              if(!success_file){

                cat('Right now ploidy equal to ploidy_ind\n')

                abline(v = avg_cov*stoped_ploidy_ind, lty=3, col="gray")
                rep_size <- genome_size - sum(as.numeric(portion_size))
                portion_size <- c(portion_size,rep_size)
                if(!ploidy_assign){
                    genome_size     = sum(histo_fit$V1 * (histo_fit$V2 / 1000000)) / (ploidy*avg_cov) # tetraploid size: 3374.929 Mb, haploid size: 843.7321
                    cat('ploidy for plot is ',stoped_ploidy_ind,'\n')
                }

                legend("topright",
                  pch    = rep(15, 3),
                  col    = c("gray", "gray", "gray","cyan", brewer_palette[2:stoped_ploidy_ind],"gray", alpha("orangered", 0.8)),
                  legend = c(paste("Raw", sep=""),
                  paste("Haploid-peak: ", avg_cov, sep=""),
                  paste("Haploid-fitting:", sep=""),
                  paste("Portion-copy",1:(stoped_ploidy_ind+1),": ",paste(round(portion_size, digits = 2), sep=","), " Mb", sep=""),
                  paste("Haploid GSE (avg.): ", round(genome_size, digits = 2), " Mb", sep="")
                  ),
                  cex    = 0.8,
                  horiz  = F,
                  box.col="NA")

                dev.off()

                next
              }
              else{
                abline(v = avg_cov*ploidy_ind, lty=3, col="gray")

                histo_raw <- histo_tmp

                ## estimate the size of haploid portion
                non_hap_size <- sum(histo_het$V1 * (histo_het$V2 / 1000000)) / (avg_cov * ploidy )
                # print(histo_het[1:100,])
                portion_size <- c(portion_size,non_hap_size)
                ## select color

                # het fitting
                lines(histo_het$V1,  histo_het$V2, col=brewer_palette[ploidy_ind], lwd=3, lty=2)
                final_ploidy <- ploidy_ind
                polygon(c(histo_het$V1, rev(histo_het$V1)),
                c(histo_het$V2, rep(1, length(histo_het$V2))),
                col = alpha(brewer_palette[ploidy_ind], 0.3), border = NA)
                # lines(histo_raw$V1,  histo_raw$V2, col="black", lwd=3, lty=2)
              }



            }
            if(ploidy_ind == ploidy){
                cat('Right now ploidy equal to ploidy_ind\n')
                rep_size <- genome_size - sum(portion_size)
                portion_size <- c(portion_size,rep_size)
                if(!ploidy_assign){
                    genome_size     = sum(histo_fit$V1 * (histo_fit$V2 / 1000000)) / (ploidy*avg_cov) # tetraploid size: 3374.929 Mb, haploid size: 843.7321
                    cat('ploidy for plot is ',final_ploidy,'\n')
                }
                legend("topright",
                  pch    = rep(15, 3),
                  col    = c("gray", "gray", "gray","cyan", brewer_palette[2:ploidy_ind],"gray", alpha("orangered", 0.8)),
                  legend = c(paste("Raw", sep=""),
                  paste("Haploid-peak: ", avg_cov, sep=""),
                  paste("Haploid-fitting:", sep=""),
                  paste("Portion-copy",1:(ploidy+1),": ",paste(round(portion_size, digits = 2), sep=","), " Mb", sep=""),
                  paste("Haploid GSE (avg.): ", round(genome_size, digits = 2), " Mb", sep="")
                  ),
                  cex    = 0.8,
                  horiz  = F,
                  box.col="NA")

                dev.off()
            }
          }





	end_time <- proc.time()
	execution_time <- end_time - start_time
	print(paste("Program running time:", sum(execution_time)))
        # png file

        png(paste(path, output_dir ,sample,"_hap_genome_size_est.png", sep=""),  width = 1200, height = 800, res = 200)
        #png(paste(path, output_dir ,sample,"_hap_genome_size_est.png", sep=""),  height=4, width=7.08661)
        par(family = "Helvetica")
        if(!success_file){
          ploidy_raw = ploidy
          ploidy = stoped_ploidy_ind
        }
        for (ploidy_ind in seq(ploidy)){
          if(ploidy_ind == 1){
            exp_hom_temp = exp_hom
            ####
            histo_raw <- read.table(paste(path, sample, sep="") ) # from jellyfish
            if((het_pos - start_pos) < 8){
                histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,"_", ploidy_ind,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R
            }
            else{
                histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R
            }
            #print(histo_raw[1:100, ])#
            #print(histo_het[1:100, ])
            # len       <- length(histo_raw$V1)
            # find the haplotype average kmer coverage.
            happeak_index <- which.max(histo_het$V2) # coverage at het-peak 14, x, ...
            #cat("hap peak at", happeak_index, '\n')
            #cat("left peak at", range_left, '\n')
            # range_avg_cov <- (happeak_index-range_left):(happeak_index+range_right)
            #avg_cov       <- round(sum(histo_raw[range_avg_cov, 1] * histo_raw[range_avg_cov, 2] / sum(histo_raw[range_avg_cov, 2])), digits = 2)

            # find the valley on the left of het region
            valley_index  <- which.min(histo_raw[1:happeak_index, 2])
            #cat(valley_index, '\n')
            #
            # make new histogram with fittings at het-kmers
            histo_fit                    <- histo_raw
            histo_fit[1:valley_index, 2] <- histo_het[1:valley_index, 2]

            if((het_pos - start_pos) < 7){
                histo_fit[1:het_pos, 2] <- histo_het[1:het_pos, 2]
            }
            #
            plot(histo_raw$V1,
                histo_raw$V2,
                xlim=c(1, xlimit),
                ylim=c(1, ylimit),
                xlab="Coverage of k-mers",
                ylab="Frequency",
                type="b",
                col ="gray",
                lwd=1,
                cex=1.5,
                frame.plot = FALSE,
                cex.axis = 0.6)
            # final used for gse
            lines(histo_fit$V1,  histo_fit$V2, col=alpha("orangered", 0.8), lwd=6)
            # het fitting
            lines(histo_het$V1,  histo_het$V2, col="cyan", lwd=3, lty=2)

            polygon(c(histo_raw$V1, rev(histo_raw$V1)),
              c(histo_raw$V2, rep(1, length(histo_raw$V2))),
              col = alpha("gray", 0.5), border = NA)

            polygon(c(histo_fit$V1, rev(histo_fit$V1)),
              c(histo_fit$V2, rep(1, length(histo_fit$V2))),
              col = alpha("orangered", 0.3), border = NA)


            polygon(c(histo_het$V1, rev(histo_het$V1)),
              c(histo_het$V2, rep(1, length(histo_het$V2))),
              col = alpha("cyan", 0.3), border = NA)
            #
            abline(v = avg_cov*1, lty=3, col="gray")

            #
            #

            #
            if(ploidy == 1){
              dev.off()
              break
            }
          }
          else{
            ## update histo_file

            histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,"_",ploidy_ind,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R
            # het fitting
            # lines(histo_het$V1,  histo_het$V2, col="cyan", lwd=3, lty=2)
            lines(histo_het$V1,  histo_het$V2, col=brewer_palette[ploidy_ind], lwd=3, lty=2)
            abline(v = avg_cov*ploidy_ind, lty=3, col="gray")
            polygon(c(histo_het$V1, rev(histo_het$V1)),
              c(histo_het$V2, rep(1, length(histo_het$V2))),
              col = alpha(brewer_palette[ploidy_ind], 0.3), border = NA)
          }
          if(ploidy_ind == ploidy){
              legend("topright",
                  pch    = rep(15, 3),
                  col    = c("gray", "gray", "gray","cyan", brewer_palette[2:ploidy_ind],"gray", alpha("orangered", 0.8)),
                  legend = c(paste("Raw", sep=""),
                  paste("Haploid-peak: ", avg_cov, sep=""),
                  paste("Haploid-fitting:", sep=""),
                  paste("Portion-copy",1:(ploidy_ind+1),": ",paste(round(portion_size, digits = 2), sep=","), " Mb", sep=""),
                  paste("Haploid GSE (avg.): ", round(genome_size, digits = 2), " Mb", sep="")
                  ),
                  cex    = 0.6,
                  horiz  = F,
                  box.col="NA")

              dev.off()
          }
        }
      }
      else{
          cat("het_pos is ",het_pos," enter het_pos > 80 process, rescaling... \n")
          if (het_pos >= 80) {
            scaled_factor <- round(het_pos / 60, digits=1)
          }
          #scaled_factor <- round(het_pos / 50)
          cat("scaled_factor is ",scaled_factor,"\n")
	  histo_raw <- read.table(paste(path, sample, sep="") ) # from jellyfish
          if((het_pos - start_pos) < 8){
                histo_org <- histo_raw
                histo_org[1:(start_pos-1),2] <- 1
                histo_data_scaled <- data.frame(
                  V1 = histo_org$V1 / scaled_factor,
                  V2 = histo_org$V2 * scaled_factor
                )
          }
          else{
              histo_data_scaled <- data.frame(
                V1 = histo_data$V1 / scaled_factor,
                V2 = histo_data$V2 * scaled_factor
              )
          }

          ## add on 20231017 Interpolation on histo_data_scaled
          interpolated_data_scaled <- approx(histo_data_scaled$V1, histo_data_scaled$V2, xout = floor(histo_data_scaled$V1))
          result_scaled <- data.frame(
            V1 = unique(interpolated_data_scaled$x),
            V2 = interpolated_data_scaled$y[!duplicated(interpolated_data_scaled$x)]
          )
          result_scaled <- na.omit(result_scaled)

          # histo_data_scaled_sel <- histo_data_scaled[histo_data_scaled[,1] %% 1 == 0,]
          # colnames(histo_data_scaled_sel) <- c("V1","V2")
          write.table(file=paste0(path,sample,"_scaled"), x=result_scaled, quote=F, row.names=F, col.names=F)
          sample <- paste0(sample,"_scaled")

          ## add main code
          portion_size <- c()
          col_list <- c()
          success_file <- TRUE
          continue_reverse_flag = FALSE
          for (ploidy_ind in seq(ploidy)){
          if(ploidy_ind == 1){
            exp_hom_temp = exp_hom
            avg_cov = 0
            cat('ploidy index : ',ploidy_ind)
            histo_raw <- read.table(paste(path, sample, sep="") ) # from jellyfish
            for(left_fit_ratio in left_fit_ratio_list){
                fit_para_save_list = findGSE_sp(histo = paste(path, sample, sep = ""), sizek = sizek , outdir = paste(path, output_dir, sep = ""), exp_hom = exp_hom_temp,ploidy_ind=ploidy_ind, avg_cov=avg_cov,left_fit_ratio=left_fit_ratio)
            #png(paste("apple_hap_genome_size_est.jpg", sep=""),  height=400, width=708)
            ####
                het_xfit_right_save <- fit_para_save_list[[1]]
                hom_peak_pos_save <- fit_para_save_list[[2]]
                het_peak_pos_save <- fit_para_save_list[[3]]
                histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R
                cat("rescaling procedure start...\n")
                histo_raw_rescale_raw <- histo_data
                #browser()
                # histo_het_rescale_tmp <- as.data.frame(cbind(histo_het$V1*scaled_factor,histo_het$V2/scaled_factor))

                # histo_het_rescale <- histo_het_rescale_tmp[histo_het_rescale_tmp[,1] %% 1 == 0,]

                # ## interpolation
                # min_v1 <- min(histo_het_rescale$V1)
                # max_v1 <- max(histo_het_rescale$V1)

                # # 创建一个包含从最小值到最大值的连续整数的数据框
                # full_v1_range <- data.frame(V1 = seq(min_v1, max_v1))

                # # 使用 approx 函数进行插值，生成新的 V2 值
                # interpolated_v2 <- approx(histo_het_rescale$V1, histo_het_rescale$V2, xout = full_v1_range$V1)$y

                # # 创建包含插值后的数据的数据框
                # histo_het_rescale_interpolated <- data.frame(V1 = full_v1_range$V1, V2 = interpolated_v2)
                ## rescale back
                histo_data_rescaled_tmp <- data.frame(
                  V1 = histo_het$V1 * scaled_factor,
                  V2 = histo_het$V2 / scaled_factor
                )
                interpolated_data_rescaled <- approx(histo_data_rescaled_tmp$V1, histo_data_rescaled_tmp$V2, xout = floor(histo_data_rescaled_tmp$V1))
                result_rescaled <- data.frame(
                  V1 = unique(interpolated_data_rescaled$x),
                  V2 = interpolated_data_rescaled$y[!duplicated(interpolated_data_rescaled$x)]
                )

                # Remove NA values from result_rescaled


                result_rescaled <- na.omit(result_rescaled)



                prepend_data <- data.frame(
                  V1 = seq(1, result_rescaled$V1[1]-1, by=1),
                  V2 = rep(0, times=length(result_rescaled$V1[1]-1))
                )

                # Combine the dataframes
                result_rescaled <- rbind(prepend_data, result_rescaled)
                ## write rescale histo file
                histo_het_rescale <- result_rescaled

                write.table(histo_het_rescale, file = paste0(path, output_dir,sample,"_",ploidy_ind,"_rescale.histo"), sep = ' ', col.names = F, row.names = F, quote = F)



                happeak_index_rescale <- which.max(histo_het_rescale$V2) # coverage at het-peak 14, x, ...
                happeak_index_rescale <- histo_het_rescale[happeak_index_rescale,1]
                range_avg_cov_rescale <- (happeak_index_rescale-range_left):(happeak_index_rescale+range_right)
                avg_cov_rescale       <- round(sum(histo_raw_rescale_raw[range_avg_cov_rescale, 1] * histo_raw_rescale_raw[range_avg_cov_rescale, 2] / sum(histo_raw_rescale_raw[range_avg_cov_rescale, 2])), digits = 2)
                cat("rescaling procedure end...\n")

                final_fit_peak_pos_ind <- which.max(histo_het_rescale$V2)
                final_fit_peak_pos <-  histo_het_rescale[final_fit_peak_pos_ind,1]
                distance_fit_two_line <- abs(avg_cov_rescale*ploidy_ind - final_fit_peak_pos)
                cat('difference of cyan peaks: ',final_fit_peak_pos,', gray peak: ',avg_cov_rescale*ploidy_ind,'\n')
                if(distance_fit_two_line <= 5){
                      cat(left_fit_ratio,' meet the requirement, continue...\n')
                      break
                }

            }
             ## re-scale back


            ## end re-scale back

            ## estimate the size of haploid portion

            #browser()
            hap_size <- sum(histo_het_rescale$V1 * histo_het_rescale$V2 / 1000000) / (avg_cov_rescale * ploidy)
            portion_size <- hap_size

            ## draw ori_plot
            len       <- length(histo_raw_rescale_raw$V1)
            # find the haplotype average kmer coverage.
            happeak_index <- which.max(histo_het$V2) # coverage at het-peak 14, x, ...
            happeak_index <- histo_het[happeak_index,1]
            cat("hap peak at", happeak_index, '\n')
            #cat("left peak at", range_left, '\n')
            if(happeak_index <= range_left){
              range_avg_cov <- 1:(happeak_index+range_right)
            }
            else{
              range_avg_cov <- (happeak_index-range_left):(happeak_index+range_right)
            }

            cat('happeak_index is: ',happeak_index,' range_left is ',range_left,'\n')
            avg_cov       <- round(sum(histo_raw[range_avg_cov, 1] * histo_raw[range_avg_cov, 2] / sum(histo_raw[range_avg_cov, 2])), digits = 2)

            # find the valley on the left of het region
            # cat('start estimate haploid portion333...: ',happeak_index_rescale,'eee ',length(histo_raw_rescale_raw[1:happeak_index_rescale, 2]),'\n')
            valley_index  <- which.min(histo_raw_rescale_raw[1:happeak_index_rescale, 2])


            #
            # make new histogram with fittings at het-kmers
            histo_fit                    <- histo_raw_rescale_raw

            ## add 0 to histo_fit
            # Find the row number in histo_het_rescale where the first column value is valley_index
            valley_index_k <- which.min(abs(histo_het_rescale[,1] - valley_index))


            cat('valley_index is: ',valley_index, ', and valley_index_k is ',valley_index_k,'\n')
            # Calculate the difference
            diff_index <- valley_index - valley_index_k

            # Construct the new vector to assign to histo_fit
            new_values <- c(rep(0, diff_index), histo_het_rescale[1:valley_index_k, 2])
            cat('histo_fit the first ',valley_index,' value is: ',new_values,'\n')
            # Assign to histo_fit
            histo_fit[1:valley_index, 2] <- new_values

            if((het_pos - start_pos) < 7){
                histo_fit[1:het_pos, 2] <- histo_het_rescale[1:het_pos, 2]
            }

            #histo_fit[1:valley_index, 2] <- histo_het_rescale[1:valley_index, 2]
            #
            if(xlimit == -1 ){
                 xlimit = avg_cov_rescale * ploidy * 2
              }
              if(ylimit == -1){
                ylimit = max(histo_het_rescale) * 1.1
              }
            plot(histo_raw_rescale_raw$V1,
                histo_raw_rescale_raw$V2,
                xlim=c(1, xlimit),
                ylim=c(1, ylimit),
                xlab="Coverage of k-mers",
                ylab="Frequency",
                type="b",
                col ="gray",
                lwd=1,
                cex=1.5,
                frame.plot = F)
            # final used for gse
            lines(histo_fit$V1,  histo_fit$V2, col=alpha("orangered", 0.8), lwd=6)

            # het fitting
            lines(histo_het_rescale$V1,  histo_het_rescale$V2, col="cyan", lwd=3, lty=2)

            polygon(c(histo_raw_rescale_raw$V1, rev(histo_raw_rescale_raw$V1)),
              c(histo_raw_rescale_raw$V2, rep(1, length(histo_raw_rescale_raw$V2))),
              col = alpha("gray", 0.5), border = NA)

            polygon(c(histo_fit$V1, rev(histo_fit$V1)),
              c(histo_fit$V2, rep(1, length(histo_fit$V2))),
              col = alpha("orangered", 0.3), border = NA)


            polygon(c(histo_het_rescale$V1, rev(histo_het_rescale$V1)),
              c(histo_het_rescale$V2, rep(1, length(histo_het_rescale$V2))),
              col = alpha("cyan", 0.3), border = NA)
            #
            abline(v = avg_cov_rescale*1, lty=3, col="gray")

            #
            genome_size_raw = sum(histo_raw_rescale_raw[valley_index:len, 1] * histo_raw_rescale_raw[valley_index:len, 2] / 1000000) / (ploidy*avg_cov_rescale)
            genome_size     = sum(histo_fit$V1 * (histo_fit$V2 / 1000000)) / (ploidy*avg_cov_rescale) # tetraploid size: 3374.929 Mb, haploid size: 843.7321
            #

            #
            if(ploidy == 1){
              dev.off()
              break
            }



          }
          else{
            if(!success_file){
                next
            }
            ## update histo_file
            histo_tmp = histo_raw
            het_length = dim(histo_het)[1]

            histo_tmp[1:het_length,2] = abs(histo_tmp[1:het_length,2] - histo_het[,2])

            ## add reverse data part


            target_start_index <- round(hom_peak_pos_save)
            target_end_index <- round(hom_peak_pos_save) * 2 - round(het_peak_pos_save)
            if(ploidy_ind > 2){
                max_het_ind <- which(histo_tmp[,2] == max(histo_tmp[round(hom_peak_pos_save):length(histo_tmp$V1),2]))
                cat('max_het_ind is : ',max_het_ind,'\n')
                if (max_het_ind && round(hom_peak_pos_save) < max_het_ind && (max_het_ind - round(hom_peak_pos_save)) > 3 ){
                    cat('match latter peak pattern, reverse data in ... ploidy ',ploidy_ind,'\n')
                    histo_tmp_raw <- histo_tmp
                    histo_tmp_raw[target_start_index:length(histo_tmp$V1),2] <- 1
                    histo_tmp_raw[target_start_index:target_end_index, 2] <- rev(histo_tmp_raw[round(het_peak_pos_save):round(hom_peak_pos_save), 2])
                    write.table(histo_tmp_raw, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = F)
                    #lines(histo_tmp_raw$V1,  histo_tmp_raw$V2, col="black", lwd=3)
                    continue_reverse_flag = TRUE
                }
                else{
                  if(continue_reverse_flag){
                  histo_tmp_raw <- histo_tmp
                  # histo_tmp_raw[which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))):length(histo_tmp_raw$V1),2] <- 1
                  histo_tmp_raw[target_start_index:length(histo_tmp$V1),2] <- 1
                  histo_tmp_raw[target_start_index:target_end_index, 2] <- rev(histo_tmp_raw[round(het_peak_pos_save):round(hom_peak_pos_save), 2])

                  #cat('pppploidy_ind range: ',which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))),'\n')
                  write.table(histo_tmp_raw, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = F)

                  }
                  else{
                  histo_tmp_raw <- histo_tmp
                  # histo_tmp_raw[which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))):length(histo_tmp_raw$V1),2] <- 1
                  #histo_tmp_raw[target_start_index:length(histo_tmp$V1),2] <- 1
                  #histo_tmp_raw[target_start_index:target_end_index, 2] <- rev(histo_tmp_raw[round(het_peak_pos_save):round(hom_peak_pos_save), 2])

                  #cat('pppploidy_ind range: ',which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))),'\n')
                  write.table(histo_tmp_raw, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = F)
                  #lines(histo_tmp_raw$V1,  histo_tmp_raw$V2, col="black", lwd=3)
                  }

                }

            }
            else{
                write.table(histo_tmp, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = F)
            }
            ## end add reverse data part

            # write.table(histo_tmp, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = F)

            exp_hom_temp = avg_cov * (ploidy_ind+1) + avg_cov * 0.1
            histo_raw <- histo_tmp
            cat('exp_hom_temp: ',exp_hom_temp,' \n')
            for(left_fit_ratio in left_fit_ratio_list){
                  cat('Now is fitting left fit ratio: ',left_fit_ratio,'\n')
                  fit_para_save_list = findGSE_sp(histo = paste(path, sample,"_",ploidy_ind, sep = ""), sizek = sizek , outdir = paste(path, output_dir, sep = ""), exp_hom = exp_hom_temp,ploidy_ind=ploidy_ind,avg_cov=avg_cov,left_fit_ratio=left_fit_ratio)
                  het_xfit_right_save <- fit_para_save_list[[1]]
                  hom_peak_pos_save <- fit_para_save_list[[2]]
                  het_peak_pos_save <- fit_para_save_list[[3]]
                  histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,"_",ploidy_ind,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R
                  # het fitting rescale
                  #lines(histo_het$V1,  histo_het$V2, col="deepskyblue", lwd=3, lty=2)
                  cat('Now scale factor is : ',scaled_factor,'\n')
                  # histo_het_rescale_tmp <- as.data.frame(cbind(histo_het$V1*scaled_factor,histo_het$V2/scaled_factor))
                  # histo_het_rescale <- histo_het_rescale_tmp[histo_het_rescale_tmp[,1] %% 1 == 0,]
                  # colnames(histo_het_rescale) <- c("V1","V2")

                  # ## interpolation
                  # min_v1 <- min(histo_het_rescale$V1)
                  # max_v1 <- max(histo_het_rescale$V1)

                  # # 创建一个包含从最小值到最大值的连续整数的数据框
                  # full_v1_range <- data.frame(V1 = seq(min_v1, max_v1))

                  # # 使用 approx 函数进行插值，生成新的 V2 值
                  # interpolated_v2 <- approx(histo_het_rescale$V1, histo_het_rescale$V2, xout = full_v1_range$V1)$y

                  # # 创建包含插值后的数据的数据框
                  # histo_het_rescale_interpolated <- data.frame(V1 = full_v1_range$V1, V2 = interpolated_v2)

                  # histo_het_rescale <- histo_het_rescale_interpolated

                  histo_data_rescaled_tmp <- data.frame(
                    V1 = histo_het$V1 * scaled_factor,
                    V2 = histo_het$V2 / scaled_factor
                  )
                  interpolated_data_rescaled <- approx(histo_data_rescaled_tmp$V1, histo_data_rescaled_tmp$V2, xout = floor(histo_data_rescaled_tmp$V1))
                  result_rescaled <- data.frame(
                    V1 = unique(interpolated_data_rescaled$x),
                    V2 = interpolated_data_rescaled$y[!duplicated(interpolated_data_rescaled$x)]
                  )

                  # Remove NA values from result_rescaled
                  result_rescaled <- na.omit(result_rescaled)

                  histo_het_rescale <- result_rescaled
                  ## write rescale histo_het_rescale
                  write.table(histo_het_rescale, file = paste0(path, output_dir,sample,"_",ploidy_ind,"_rescale.histo"), sep = ' ', col.names = F, row.names = F, quote = F)
                  all_slope_data <- abs(diff(diff(histo_het$V2))[(round(avg_cov)*ploidy_ind-10):(round(avg_cov)*ploidy_ind+10)])
                  if(!all(all_slope_data > 100)){
                    #lines(histo_het_rescale$V1,  histo_het_rescale$V2, col="deepskyblue", lwd=3, lty=2)

                    cat('index is :',ploidy_ind,' slope is ',all_slope_data)
                    success_file <- FALSE
                    stoped_ploidy_ind = ploidy_ind - 1
                    break
                  }
                  final_fit_peak_pos_ind <- which.max(histo_het_rescale$V2)

                  final_fit_peak_pos <-  histo_het_rescale[final_fit_peak_pos_ind,1]
                  distance_fit_two_line <- abs(avg_cov_rescale*ploidy_ind - final_fit_peak_pos)
                  cat('difference of cyan peaks: ',final_fit_peak_pos,', gray peak: ',avg_cov_rescale*ploidy_ind,'\n')
                  if(distance_fit_two_line <= 5){
                      cat(left_fit_ratio,' meet the requirement, continue...\n')
                      break
                  }
                  cat('Do Not meet requirement, try another parameter...\n')
            }

            if(!success_file){

                cat('Right now ploidy equal to ploidy_ind\n')

                abline(v = avg_cov_rescale*stoped_ploidy_ind, lty=3, col="gray")
                rep_size <- genome_size - sum(portion_size)
                portion_size <- c(portion_size,rep_size)
                legend("topright",
                  pch    = rep(15, 3),
                  col    = c("gray", "gray", "gray","cyan", brewer_palette[2:stoped_ploidy_ind],"gray", alpha("orangered", 0.8)),
                  legend = c(paste("Raw", sep=""),
                  paste("Haploid-peak: ", avg_cov_rescale, sep=""),
                  paste("Haploid-fitting:", sep=""),
                  paste("Portion-copy",1:(stoped_ploidy_ind+1),": ",paste(round(portion_size, digits = 2), sep=","), " Mb", sep=""),
                  paste("Haploid GSE (avg.): ", round(genome_size, digits = 2), " Mb", sep="")
                  ),
                  cex    = 0.8,
                  horiz  = F,
                  box.col="NA")

                dev.off()

                next
              }
              else{
                non_hap_size <- sum(histo_het_rescale$V1 * histo_het_rescale$V2 / 1000000) / (avg_cov_rescale * ploidy )

                portion_size <- c(portion_size,non_hap_size)
                abline(v = avg_cov_rescale*ploidy_ind, lty=3, col="gray")
                # end het fitting rescale
                lines(histo_het_rescale$V1,  histo_het_rescale$V2, col=brewer_palette[ploidy_ind], lwd=3, lty=2)
                #lines(histo_raw$V1,  histo_raw$V2, col="black", lwd=3, lty=2)
                polygon(c(histo_het_rescale$V1, rev(histo_het_rescale$V1)),
                  c(histo_het_rescale$V2, rep(1, length(histo_het_rescale$V2))),
                  col = alpha(brewer_palette[ploidy_ind], 0.3), border = NA)
              }
             ## estimate the size of haploid portion

          }
          if(ploidy_ind == ploidy){
              # rep_start_ind <- round(avg_cov_rescale * (ploidy + 0.5))
              #rep_size <- sum(histo_fit[rep_start_ind:length(histo_fit$V1),1] * histo_fit[rep_start_ind:length(histo_fit$V1),2] / 1000000) / (ploidy*avg_cov_rescale)
              rep_size <- genome_size - sum(portion_size)
              portion_size <- c(portion_size,rep_size)
              legend("topright",
                  pch    = rep(15, 3),
                  col    = c("gray", "gray", "gray","cyan", brewer_palette[2:ploidy_ind],"gray", alpha("orangered", 0.8)),
                  legend = c(paste("Raw", sep=""),
                  paste("Haploid-peak: ", avg_cov_rescale, sep=""),
                  paste("Haploid-fitting", sep=""),
                  paste("Portion-copy",1:(ploidy_ind+1),": ",paste(round(portion_size, digits = 2), sep=","), " Mb", sep=""),
                  paste("Haploid GSE (avg.): ", round(genome_size, digits = 2), " Mb", sep="")
                  ),
                  cex    = 0.8,
                  horiz  = F,
                  box.col="NA")
              dev.off()
          }
        }

      ## end main code

      end_time <- proc.time()
      execution_time <- end_time - start_time
      print(paste("Program running time:", sum(execution_time)))
      ## png file draw


      png(paste(path, output_dir ,sub("_scaled$", "", sample),"_hap_genome_size_est.png", sep=""),  width = 1200, height = 800, res = 200)
      par(family = "Helvetica")
      if(!success_file){
          ploidy_raw = ploidy
          ploidy = stoped_ploidy_ind
      }
      for (ploidy_ind in seq(ploidy)){
          if(ploidy_ind == 1){

            avg_cov = 0
            cat('ploidy index : ',ploidy_ind)
            #png(paste("apple_hap_genome_size_est.jpg", sep=""),  height=400, width=708)
            ####
            histo_raw <- read.table(paste(path, sub("_scaled$", "", sample), sep="") ) # from jellyfish
            histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R
            ## re-scale back
            cat("rescaling procedure start...\n")
            # histo_raw_rescale_raw <- histo_data
            # #browser()
            # histo_het_rescale_tmp <- as.data.frame(cbind(histo_het$V1*scaled_factor,histo_het$V2/scaled_factor))

            # histo_het_rescale <- histo_het_rescale_tmp[histo_het_rescale_tmp[,1] %% 1 == 0,]

            # ## interpolation
            # min_v1 <- min(histo_het_rescale$V1)
            # max_v1 <- max(histo_het_rescale$V1)

            # # 创建一个包含从最小值到最大值的连续整数的数据框
            # full_v1_range <- data.frame(V1 = seq(min_v1, max_v1))

            # # 使用 approx 函数进行插值，生成新的 V2 值
            # interpolated_v2 <- approx(histo_het_rescale$V1, histo_het_rescale$V2, xout = full_v1_range$V1)$y

            # # 创建包含插值后的数据的数据框
            # histo_het_rescale_interpolated <- data.frame(V1 = full_v1_range$V1, V2 = interpolated_v2)


            # histo_het_rescale <- histo_het_rescale_interpolated
            histo_data_rescaled_tmp <- data.frame(
                    V1 = histo_het$V1 * scaled_factor,
                    V2 = histo_het$V2 / scaled_factor
                  )
            interpolated_data_rescaled <- approx(histo_data_rescaled_tmp$V1, histo_data_rescaled_tmp$V2, xout = floor(histo_data_rescaled_tmp$V1))
            result_rescaled <- data.frame(
                    V1 = unique(interpolated_data_rescaled$x),
                    V2 = interpolated_data_rescaled$y[!duplicated(interpolated_data_rescaled$x)]
                  )

            # Remove NA values from result_rescaled
            result_rescaled <- na.omit(result_rescaled)

            histo_het_rescale <- result_rescaled



            happeak_index_rescale <- which.max(histo_het_rescale$V2) # coverage at het-peak 14, x, ...
            range_avg_cov_rescale <- (happeak_index_rescale-range_left):(happeak_index_rescale+range_right)
            #avg_cov_rescale       <- round(sum(histo_raw_rescale_raw[range_avg_cov_rescale, 1] * histo_raw_rescale_raw[range_avg_cov_rescale, 2] / sum(histo_raw_rescale_raw[range_avg_cov_rescale, 2])), digits = 2)
            cat("rescaling procedure end...\n")

            ## end re-scale back




            # find the haplotype average kmer coverage.
            happeak_index <- which.max(histo_het$V2) # coverage at het-peak 14, x, ...
            #cat("hap peak at", happeak_index, '\n')
            #cat("left peak at", range_left, '\n')
            range_avg_cov <- (happeak_index-range_left):(happeak_index+range_right)
            avg_cov       <- round(sum(histo_raw[range_avg_cov, 1] * histo_raw[range_avg_cov, 2] / sum(histo_raw[range_avg_cov, 2])), digits = 2)

            # find the valley on the left of het region
            valley_index  <- which.min(histo_raw_rescale_raw[1:happeak_index_rescale, 2])

            #cat(valley_index, '\n')
            #


            #
            plot(histo_raw_rescale_raw$V1,
                histo_raw_rescale_raw$V2,
                xlim=c(1, xlimit),
                ylim=c(1, ylimit),
                xlab="Coverage of k-mers",
                ylab="Frequency",
                type="b",
                col ="gray",
                lwd=1,
                cex=1.5,
                frame.plot = FALSE,
                cex.axis = 0.6)
            # final used for gse
            lines(histo_fit$V1,  histo_fit$V2, col=alpha("orangered", 0.8), lwd=6)
            # het fitting
            lines(histo_het_rescale$V1,  histo_het_rescale$V2, col="cyan", lwd=3, lty=2)

            polygon(c(histo_raw_rescale_raw$V1, rev(histo_raw_rescale_raw$V1)),
              c(histo_raw_rescale_raw$V2, rep(1, length(histo_raw_rescale_raw$V2))),
              col = alpha("gray", 0.5), border = NA)

            polygon(c(histo_fit$V1, rev(histo_fit$V1)),
              c(histo_fit$V2, rep(1, length(histo_fit$V2))),
              col = alpha("orangered", 0.3), border = NA)


            polygon(c(histo_het_rescale$V1, rev(histo_het_rescale$V1)),
              c(histo_het_rescale$V2, rep(1, length(histo_het_rescale$V2))),
              col = alpha("cyan", 0.3), border = NA)
            #
            abline(v = avg_cov_rescale*1, lty=3, col="gray")

            #


            #
            if(ploidy == 1){
              dev.off()
              break
            }



          }
          else{
            ## update histo_file

            abline(v = avg_cov_rescale*ploidy_ind, lty=3, col="gray")
            histo_het_rescale <- read.table(paste0(path, output_dir,sample,"_",ploidy_ind,"_rescale.histo"), sep = ' ')
            # end het fitting rescale
            lines(histo_het_rescale$V1,  histo_het_rescale$V2, col=brewer_palette[ploidy_ind], lwd=3, lty=2)
            polygon(c(histo_het_rescale$V1, rev(histo_het_rescale$V1)),
              c(histo_het_rescale$V2, rep(1, length(histo_het_rescale$V2))),
              col = alpha(brewer_palette[ploidy_ind], 0.3), border = NA)
            #lines(histo_raw$V1,  histo_raw$V2, col="black", lwd=3, lty=2)
          }
          if(ploidy_ind == ploidy){
              # rep_start_ind <- round(avg_cov_rescale * (ploidy + 0.5))
              #rep_size <- sum(histo_fit[rep_start_ind:length(histo_fit$V1),1] * histo_fit[rep_start_ind:length(histo_fit$V1),2] / 1000000) / (ploidy*avg_cov_rescale)

              legend("topright",
                  pch    = rep(15, 3),
                  col    = c("gray", "gray", "gray","cyan", brewer_palette[2:ploidy_ind],"gray", alpha("orangered", 0.8)),
                  legend = c(paste("Raw", sep=""),
                  paste("Haploid-peak: ", avg_cov_rescale, sep=""),
                  paste("Haploid-fitting", sep=""),
                  paste("Portion-copy",1:(ploidy_ind+1),": ",paste(round(portion_size, digits = 2), sep=","), " Mb", sep=""),
                  paste("Haploid GSE (avg.): ", round(genome_size, digits = 2), " Mb", sep="")
                  ),
                  cex    = 0.6,
                  horiz  = F,
                  box.col="NA")
              dev.off()
          }
        }
	}
  }






  #    # png file

  #   png(paste(path, output_dir ,sample,"_hap_genome_size_est.png", sep=""),  width = 800, height = 600)

  #   for (ploidy_ind in seq(ploidy)){
  #     if(ploidy_ind == 1){
  #       exp_hom_temp = exp_hom
  #       ####
  #       histo_raw <- read.table(paste(path, sample, sep="") ) # from jellyfish
  #       histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R
  #       #print(histo_raw[1:100, ])#
  #       #print(histo_het[1:100, ])
  #       len       <- length(histo_raw$V1)
  #       # find the haplotype average kmer coverage.
  #       happeak_index <- which.max(histo_het$V2) # coverage at het-peak 14, x, ...
  #       #cat("hap peak at", happeak_index, '\n')
  #       #cat("left peak at", range_left, '\n')
  #       range_avg_cov <- (happeak_index-range_left):(happeak_index+range_right)
  #       avg_cov       <- round(sum(histo_raw[range_avg_cov, 1] * histo_raw[range_avg_cov, 2] / sum(histo_raw[range_avg_cov, 2])), digits = 2)
  #       # find the valley on the left of het region
  #       valley_index  <- which.min(histo_raw[1:happeak_index, 2])
  #       #cat(valley_index, '\n')
  #       #
  #       # make new histogram with fittings at het-kmers
  #       histo_fit                    <- histo_raw
  #       histo_fit[1:valley_index, 2] <- histo_het[1:valley_index, 2]
  #       #
  #       plot(histo_raw$V1,
  #           histo_raw$V2,
  #           xlim=c(1, xlimit),
  #           ylim=c(1, ylimit),
  #           xlab="Coverage of k-mers",
  #           ylab="Frequency",
  #           type="b",
  #           col ="gray",
  #           lwd=1,
  #           cex=1.5,
  #           frame.plot = F)
  #       # final used for gse
  #       lines(histo_fit$V1,  histo_fit$V2, col=alpha("orangered", 0.8), lwd=6)
  #       # het fitting
  #       lines(histo_het$V1,  histo_het$V2, col="cyan", lwd=3, lty=2)
  #       #
  #       abline(v = avg_cov*1, lty=3, col="gray")
  #       abline(v = avg_cov*2, lty=3, col="gray")
  #       abline(v = avg_cov*3, lty=3, col="gray")
  #       abline(v = avg_cov*4, lty=3, col="gray")
  #       #
  #       #
  #       legend("topright",
  #             pch    = rep(15, 3),
  #             col    = c("gray", "gray", "cyan", alpha("orangered", 0.8)),
  #             legend = c(paste("Raw", sep=""),
  #             paste("Haploid-peak: ", avg_cov, sep=""),
  #             paste("Haploid-fitting", sep=""),
  #             paste("Haploid GSE (avg.): ", round(genome_size, digits = 2), " Mb", sep="")
  #             ),
  #             cex    = 0.8,
  #             horiz  = F,
  #             box.col="NA")
  #       #
  #       if(ploidy == 1){
  #         dev.off()
  #         break
  #       }
  #     }
  #     else{
  #       ## update histo_file
  #       histo_tmp = histo_raw
  #       het_length = dim(histo_het)[1]
  #       histo_tmp[1:het_length,2] = abs(histo_tmp[1:het_length,2] - histo_het[,2])
	# #histo_tmp[1:(ploidy_ind/(ploidy_ind+1)*avg_cov*ploidy_ind)*0.6,2] = 0
  #       write.table(histo_tmp, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = F)

  #       exp_hom_temp = avg_cov * (ploidy_ind+1) + 50
  #       findGSE_sp(histo = paste(path, sample,"_",ploidy_ind, sep = ""), sizek = sizek , outdir = paste(path, output_dir, sep = ""), exp_hom = exp_hom_temp,ploidy_ind=ploidy_ind,avg_cov=avg_cov)

  #       histo_raw <- histo_tmp
  #       histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,"_",ploidy_ind,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R
  #       # het fitting
  #       lines(histo_het$V1,  histo_het$V2, col="cyan", lwd=3, lty=2)
  #     }
  #     if(ploidy_ind == ploidy){
  #         dev.off()
  #     }
  #   }
    haploid_size_data <- data.frame()
    if (ploidy_raw > 0){
        ploidy = ploidy_raw
    }
    if (ploidy == 1) {
      haploid_size_data <- data.frame(monoploid_size = paste0(as.integer(genome_size * ploidy), " Mb"), haploid_size = paste0(as.integer(genome_size), " Mb"))
    } else if (ploidy == 2) {
      haploid_size_data <- data.frame(diploid_size = paste0(as.integer(genome_size * ploidy), " Mb"), haploid_size = paste0(as.integer(genome_size), " Mb"))
    } else if (ploidy == 3) {
      haploid_size_data <- data.frame(triploid_size = paste0(as.integer(genome_size * ploidy), " Mb"), haploid_size = paste0(as.integer(genome_size), " Mb"))
    } else if (ploidy == 4) {
      haploid_size_data <- data.frame(tetraploid_size = paste0(as.integer(genome_size * ploidy), " Mb"), haploid_size = paste0(as.integer(genome_size), " Mb"))
    } else if (ploidy == 5) {
      haploid_size_data <- data.frame(pentaploid_size = paste0(as.integer(genome_size * ploidy), " Mb"), haploid_size = paste0(as.integer(genome_size), " Mb"))
    } else if (ploidy == 6) {
      haploid_size_data <- data.frame(hexaploid_size = paste0(as.integer(genome_size * ploidy), " Mb"), haploid_size = paste0(as.integer(genome_size), " Mb"))
    } else if (ploidy == 7) {
      haploid_size_data <- data.frame(heptaploid_size = paste0(as.integer(genome_size * ploidy), " Mb"), haploid_size = paste0(as.integer(genome_size), " Mb"))
    } else if (ploidy == 8) {
      haploid_size_data <- data.frame(octoploid_size = paste0(as.integer(genome_size * ploidy), " Mb"), haploid_size = paste0(as.integer(genome_size), " Mb"))
    } else if (ploidy == 9) {
      haploid_size_data <- data.frame(nonaploid_size = paste0(as.integer(genome_size * ploidy), " Mb"), haploid_size = paste0(as.integer(genome_size), " Mb"))
    } else if (ploidy >= 10) {
      haploid_size_data <- data.frame(decaploid_size = paste0(as.integer(genome_size * ploidy), " Mb"), haploid_size = paste0(as.integer(genome_size), " Mb"))
    }
    write.csv(haploid_size_data, file = paste(path, output_dir ,sample,"_haploid_size.csv", sep=""), sep = ",",row.names=FALSE, quote=FALSE)

   }
}
