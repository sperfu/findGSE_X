# findGSEP R Package - Estimate Genome Size
# Version 1.2
#
# Description:
# The findGSEP R package provides a method for estimating genome size by fitting
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
# devtools::install_github("sperfu/findGSEP")
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
# Call the findGSEP function with specified parameters:
# findGSEP(path, samples, sizek, exp_hom, ploidy, range_left, range_right, xlimit, ylimit, output_dir)
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
source('R/utils.R')
library(scales)
library(dplyr)
# library(RColorBrewer)
library(ggplot2)
library(grDevices)

################################################## main #############################################################
#' @title Estimate genome size of polyploid species using k-mer frequencies.
#'
#' @description findGSEP is a function for multiple polyploidy
#' genome size estimation by fitting k-mer frequencies iteratively
#' with a normal distribution model.
#'
#' @description To use findGSEP, one needs to prepare a histo file,
#' which contains two tab-separated columns.
#' The first column gives frequencies at which k-mers occur in reads,
#' while the second column gives counts of such distinct k-mers.
#' Parameters k and related histo file are required for any estimation.
#'
#' @description Dependencies (R library) required: pracma, fGarch, etc. - see DESCRIPTION for details.
#'
#' @import dnorm,pdf,brewer.pal
#' @param path is the histo file location (mandatory).
#' @param samples is the histo file name (mandatory)
#' @param sizek is the size of k used to generate the histo file (mandatory).
#' K is involved in calculating heterzygosity if the genome is heterozygous.
#' @param exp_hom a rough average k-mer coverage for finding the homozygous regions.
#' In general, one can get peaks in the k-mer frequencies file, but has to
#' determine which one is for the homozygous regions, and which one is for the
#' heterozygous regions. It is optional, however, it must be provided
#' if one wants to estimate size for a heterozygous genome.
#' VALUE for exp_hom must satisfy fp < VALUE < 2*fp, where fp is the freq for homozygous peak.
#' If not provided, 0 by default assumes the genome is homozygous.
#' @param ploidy is the number of ploidy. (mandatory).
#' @param range_left is the left range for estimation, default is exp_hom*0.2, normally do not need
#' to change this. (optional).
#' @param range_right is the right range for estimation, default is exp_hom*0.2, normally do not need
#' to change this. (optional).
#' @param xlimit is the x-axis range, if not given, then it will automatically calculate a proper range,
#' normally do not need to change this. (optional).
#' @param ylimit is the y-axis range, if not given, then it will automatically calculate a proper range,
#' normally do not need to change this. (optional).
#' @param output_dir is the path to write output files (optional).
#' If not provided, by default results will be written in the folder
#' where the histo file is.
#' @return No return value, called for side effects. The function generates PDF, PNG, and CSV files in the specified output directory.
#'
#' @examples
#' \donttest{
#' test_histo <- system.file("extdata","example.histo",package = "findGSEP")
#' path <- dirname(test_histo)
#' samples <- basename(test_histo)
#' sizek <- 21
#' exp_hom <- 200
#' ploidy <- 3
#' range_left <- exp_hom*0.2
#' range_right <- exp_hom*0.2
#' xlimit <- -1
#' ylimit <- -1
#' output_dir <- ""
#'
#' findGSEP(path, samples, sizek, exp_hom, ploidy, range_left, range_right, xlimit, ylimit, output_dir)
#' }
#' @export
#'
findGSEP <- function(path, samples, sizek, exp_hom, ploidy, range_left, range_right, xlimit, ylimit ,output_dir="outfile"){
  oldpar <- par(no.readonly = TRUE) #
  on.exit(par(oldpar)) #
  start_time <- initialize_start_time()
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

  #message(range_left,'\n')
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
  meanfit_old <- 0
  sdfit_old <- 0
  for (sample in samples) {

    pdf(paste(path, output_dir ,sample,"_hap_genome_size_est.pdf", sep=""), family="Helvetica", height=4, width=7.08661) # need to update

    if(ploidy <= 2){
      message('Ploidy less than 2 ,starting...\n')
      brewer_palette <- RColorBrewer::brewer.pal(ploidy, "Set1")
      findGSE_raw(histo=paste0(path ,sample),sizek=sizek, outdir=paste0(path, output_dir), exp_hom=exp_hom)

      message('Start findGSE plot drawing, please wait...\n')
      ## read data

      histo_raw <- read.table(paste0(path ,sample))
      message('histo_raw read done...\n')
      histo_fit <- read.table(paste(path, output_dir,'v1.95.est.',
                                    sample, ".genome.size.estimated.k",
                                    min(sizek), 'to', max(sizek),".fitted_fullfit_count.txt",sep=""))

      message('histo_fit read done...\n')
      if(exp_hom != 0){
        histo_het <- read.table(paste(path, output_dir,'v1.95.est.',
                                      sample, ".genome.size.estimated.k",
                                      min(sizek), 'to', max(sizek),".fitted_hetfit_count.txt",sep=""))
        message('histo_het read done...\n')
      }
      else{
        histo_het <- read.table(paste(path, output_dir,'v1.95.est.',
                                      sample, ".genome.size.estimated.k",
                                      min(sizek), 'to', max(sizek),".fitted_hetfit_count.txt",sep=""))
        histo_first_fit_hom <-   histo_het
        # histo_first_fit_hom <- read.table(paste(path, output_dir,'v1.95.est.',
        #                   sample, ".genome.size.estimated.k",
        #                   min(sizek), 'to', max(sizek),".fitted_first_fit_hom_count.txt",sep=""))
        message('first_fit_hom read done...\n')
      }


      #load data
      data_to_save <- NULL
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
           frame.plot = FALSE)
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
      message('plot het line done...\n')

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

      end_time <- initialize_start_time()
      execution_time <- end_time - start_time
      message(paste("Program running time:", sum(execution_time)))
      ## draw grDevices::png file
      grDevices::png(paste(path, output_dir ,sample,"_hap_genome_size_est.png", sep=""),  width = 1200, height = 800, res = 200)
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
      message("ploidy id larger that 2, enter next process...",path,sample,"\n")
      histo_data <- read.table(paste(path, sample, sep = ""), header = FALSE)
      het_pos_list <- get_het_pos(histo_data)
      start_pos <- het_pos_list[[1]]
      het_pos <- het_pos_list[[2]]
      # left_fit_ratio_list <- c(0.835, 0.855, 0.825)
      left_fit_ratio_list <- c(0.835)
      ## define selected color
      # select color
      brewer_palette <- RColorBrewer::brewer.pal(ploidy, "Set1")
      if(het_pos < 80){
        message("het_pos is ",het_pos," enter het_pos < 80 process... \n")
        portion_size <- c()
        success_file <- TRUE
        continue_reverse_flag = FALSE
        scale_flag <- FALSE
        for (ploidy_ind in seq(ploidy)){
          if(ploidy_ind == 1){
            exp_hom_temp = exp_hom
            avg_cov = 0
            message('ploidy index : ',ploidy_ind)
            histo_raw <- read.table(paste(path, sample, sep="") ) # from jellyfish
            ## 1204 add cutting
            if((het_pos - start_pos) < 8){
              histo_org <- histo_raw
              histo_org[1:(start_pos-1),2] <- 1
              # pivot_value <- histo_org$V2[start_pos]
              # histo_org$V2[1:(start_pos-1)] <- pivot_value - (histo_org$V2[1:(start_pos-1)] - pivot_value)


              # histo_org[(het_pos-29+het_pos):het_pos, 2] <- rev(histo_org[round(het_pos):29, 2])

              # histo_org[30:2000, 2] <- 1
              write.table(histo_org, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
            }
            ##end 1204
            for(left_fit_ratio in left_fit_ratio_list){
              if((het_pos - start_pos) < 8){
                message("read cutting file...\n")
                fit_para_save_list = findGSE_sp(histo = paste(path, sample,"_",ploidy_ind, sep = ""), sizek = sizek , outdir = paste(path, output_dir, sep = ""), exp_hom = exp_hom_temp,ploidy_ind=ploidy_ind, avg_cov=avg_cov, left_fit_ratio=left_fit_ratio, meanfit_old=meanfit_old, sdfit_old = sdfit_old, scale_flag = scale_flag)
                het_xfit_right_save <- fit_para_save_list[[1]]
                hom_peak_pos_save <- fit_para_save_list[[2]]
                het_peak_pos_save <- fit_para_save_list[[3]]
                meanfit_old <- fit_para_save_list[[4]]
                sdfit_old <- fit_para_save_list[[5]]
                histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample, "_", ploidy_ind, ".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R
              }
              else{
                fit_para_save_list = findGSE_sp(histo = paste(path, sample, sep = ""), sizek = sizek , outdir = paste(path, output_dir, sep = ""), exp_hom = exp_hom_temp,ploidy_ind=ploidy_ind, avg_cov=avg_cov, left_fit_ratio=left_fit_ratio, meanfit_old=meanfit_old, sdfit_old = sdfit_old, scale_flag = scale_flag)
                het_xfit_right_save <- fit_para_save_list[[1]]
                hom_peak_pos_save <- fit_para_save_list[[2]]
                het_peak_pos_save <- fit_para_save_list[[3]]
                meanfit_old <- fit_para_save_list[[4]]
                sdfit_old <- fit_para_save_list[[5]]
                message('meanfit_old in ploidy 1 is:', meanfit_old,' \n')
                histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R

              }                  #grDevices::png(paste("apple_hap_genome_size_est.jpg", sep=""),  height=400, width=708)
              ####
              #print(histo_raw[1:100, ])#
              #print(histo_het[1:100, ])
              len       <- length(histo_raw$V1)
              # find the haplotype average kmer coverage.
              happeak_index <- which.max(histo_het$V2) # coverage at het-peak 14, x, ...
              #message("hap peak at", happeak_index, '\n')
              #message("left peak at", range_left, '\n')
              range_avg_cov <- (happeak_index-range_left):(happeak_index+range_right)
              avg_cov       <- round(sum(histo_raw[range_avg_cov, 1] * histo_raw[range_avg_cov, 2] / sum(histo_raw[range_avg_cov, 2])), digits = 2)
              message('avg_cov: ',avg_cov,'\n')
              final_fit_peak_pos_ind <- which.max(histo_het$V2)
              final_fit_peak_pos <-  histo_het[final_fit_peak_pos_ind,1]
              distance_fit_two_line <- abs(avg_cov*ploidy_ind - final_fit_peak_pos)
              message('difference of cyan peaks: ',final_fit_peak_pos,', gray peak: ',avg_cov*ploidy_ind,'\n')
              if(distance_fit_two_line <= 5){
                message(left_fit_ratio,' meet the requirement, continue...\n')
                break
              }
            }



            ## estimate the size of haploid portion
            hap_size <- sum(as.numeric(histo_het$V1) * as.numeric(histo_het$V2) / 1000000) / (avg_cov * ploidy)
            message('portion_size for ploidy_ind 1: ',hap_size,'\n')
            portion_size <- hap_size

            # find the valley on the left of het region
            valley_index  <- which.min(histo_raw[1:happeak_index, 2])
            #message(valley_index, '\n')
            #
            # make new histogram with fittings at het-kmers
            histo_fit                    <- histo_raw
            histo_fit[1:valley_index, 2] <- histo_het[1:valley_index, 2]
            #
            if((het_pos - start_pos) < 7){
              message('het_pos - start_pos) < 7 ohhhhh myyyy!!\n')
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
                 frame.plot = FALSE)
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
            message('difference of cyan peaks: ',final_fit_peak_pos,', gray peak: ',avg_cov,'\n')

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
            message("ploidy_ind is: ",ploidy_ind)
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
              message('max_het_ind is : ',max_het_ind,'\n')
              if (max_het_ind && round(hom_peak_pos_save) < max_het_ind && (max_het_ind - round(hom_peak_pos_save)) > 3 ){
                message('match latter peak pattern, reverse data in ... ploidy ',ploidy_ind,'\n')
                histo_tmp_raw <- histo_tmp
                histo_tmp_raw[target_start_index:length(histo_tmp$V1),2] <- 1
                histo_tmp_raw[target_start_index:target_end_index, 2] <- rev(histo_tmp_raw[round(het_peak_pos_save):round(hom_peak_pos_save), 2])
                histo_tmp_raw[target_start_index,2] <- histo_tmp_raw[target_start_index+1,2]
                write.table(histo_tmp_raw, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
                #lines(histo_tmp_raw$V1,  histo_tmp_raw$V2, col="black", lwd=3)
                continue_reverse_flag = TRUE
              }
              else{
                if(continue_reverse_flag){
                  message('continue_reverse_flag is TRUE... \n')
                  histo_tmp_raw <- histo_tmp
                  # histo_tmp_raw[which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))):length(histo_tmp_raw$V1),2] <- 1
                  # histo_tmp_raw[target_start_index:length(histo_tmp$V1),2] <- 1
                  # histo_tmp_raw[1:round(meanfit_old*(ploidy_ind-1)),2] <- 1

                  #histo_tmp_raw[target_start_index:target_end_index, 2] <- rev(histo_tmp_raw[round(het_peak_pos_save):round(hom_peak_pos_save), 2])
                  #histo_tmp_raw[target_start_index,2] <- histo_tmp_raw[target_start_index+1,2]
                  #lines(histo_tmp_raw$V1,  histo_tmp_raw$V2, col="black", lwd=3, lty=2)
                  #message('pppploidy_ind range: ',which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))),'\n')
                  write.table(histo_tmp_raw, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)

                }
                else{
                  histo_tmp_raw <- histo_tmp
                  message('This difference is ',histo_tmp[avg_cov*ploidy_ind,2] / histo_tmp[het_xfit_right_save-1,2],'\n')
                  message('The two peak difference is ',histo_tmp[avg_cov*(ploidy_ind+1),2] / histo_tmp[avg_cov*ploidy_ind,2],'\n')
                  message('This peak difference is ',histo_tmp[avg_cov*(ploidy_ind+1),2] / histo_tmp[het_xfit_right_save-1,2],'\n')
                  if(histo_tmp_raw[avg_cov*ploidy_ind,2] / histo_tmp_raw[het_xfit_right_save-1,2] < 3 && histo_tmp_raw[avg_cov*ploidy_ind,2] / histo_tmp_raw[het_xfit_right_save-1,2] > 1.6){
                    message('Difference too close, lower the histo_tmp file for ',ploidy_ind,'\n')
                    histo_tmp_raw[1:het_xfit_right_save,2] <- 1
                    #histo_tmp[avg_cov*ploidy_ind:het_xfit_right_save,2] <- round(histo_tmp[avg_cov*ploidy_ind,2])
                  }
                  # histo_tmp_raw[which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))):length(histo_tmp_raw$V1),2] <- 1
                  #histo_tmp_raw[target_start_index:length(histo_tmp$V1),2] <- 1
                  #histo_tmp_raw[target_start_index:target_end_index, 2] <- rev(histo_tmp_raw[round(het_peak_pos_save):round(hom_peak_pos_save), 2])

                  #message('pppploidy_ind range: ',which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))),'\n')
                  write.table(histo_tmp_raw, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
                  # lines(histo_tmp_raw$V1,  histo_tmp_raw$V2, col="black", lwd=3)
                }

              }

            }
            else{
              message('This difference is ',histo_tmp[avg_cov*ploidy_ind,2] / histo_tmp[het_xfit_right_save-1,2],'\n')
              message('The two peak difference is ',histo_tmp[avg_cov*(ploidy_ind+1),2] / histo_tmp[avg_cov*ploidy_ind,2],'\n')
              message('This peak difference is ',histo_tmp[avg_cov*(ploidy_ind+1),2] / histo_tmp[het_xfit_right_save-1,2],'\n')
              if(histo_tmp[avg_cov*ploidy_ind,2] / histo_tmp[het_xfit_right_save-1,2] < 3 && histo_tmp[avg_cov*ploidy_ind,2] / histo_tmp[het_xfit_right_save-1,2] > 1.6){
                message('Difference too close, lower the histo_tmp file for ',ploidy_ind,'\n')
                histo_tmp[1:het_xfit_right_save,2] <- 1
                message('histo_tmp CHANGE AGAIN!!!!\n')
                #histo_tmp[avg_cov*ploidy_ind:het_xfit_right_save,2] <- round(histo_tmp[avg_cov*ploidy_ind,2])
              }
              # lines(histo_tmp$V1,  histo_tmp$V2, col="black", lwd=3)
              write.table(histo_tmp, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
            }
            ## end add reverse data part





            ## may delete
            #write.table(histo_tmp, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = FALSE)

            exp_hom_temp = avg_cov * (ploidy_ind+1) + avg_cov * 0.1
            for(left_fit_ratio in left_fit_ratio_list){


              fit_para_save_list = findGSE_sp(histo = paste(path, sample,"_",ploidy_ind, sep = ""), sizek = sizek , outdir = paste(path, output_dir, sep = ""), exp_hom = exp_hom_temp,ploidy_ind=ploidy_ind,avg_cov=avg_cov,left_fit_ratio=left_fit_ratio, meanfit_old=meanfit_old, sdfit_old = sdfit_old,scale_flag = scale_flag)
              het_xfit_right_save <- fit_para_save_list[[1]]
              hom_peak_pos_save <- fit_para_save_list[[2]]
              het_peak_pos_save <- fit_para_save_list[[3]]
              histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,"_",ploidy_ind,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R

              all_slope_data <- abs(diff(diff(histo_het$V2))[(round(avg_cov)*ploidy_ind-10):(round(avg_cov)*ploidy_ind+10)])
              if(!all(all_slope_data > 0)){
                message('index is :',ploidy_ind,' slope is :',all_slope_data,'\n')
                success_file <- FALSE
                stoped_ploidy_ind = ploidy_ind - 1
                break
              }

              final_fit_peak_pos_ind <- which.max(histo_het$V2)
              final_fit_peak_pos <-  histo_het[final_fit_peak_pos_ind,1]
              distance_fit_two_line <- abs(avg_cov*ploidy_ind - final_fit_peak_pos)
              message('difference of cyan peaks: ',final_fit_peak_pos,', gray peak: ',avg_cov*ploidy_ind,'\n')
              if(distance_fit_two_line <= 5){
                message(left_fit_ratio,' meet the requirement, continue...\n')
                break
              }

            }
            if(!success_file){

              message('Right now ploidy equal to ploidy_ind\n')



              #abline(v = avg_cov*stoped_ploidy_ind, lty=3, col="gray")
              abline(v = meanfit_old*stoped_ploidy_ind, lty=3, col="gray")

              rep_size <- genome_size - sum(as.numeric(portion_size))
              portion_size <- c(portion_size,rep_size)
              if(!ploidy_assign){
                genome_size     = sum(histo_fit$V1 * (histo_fit$V2 / 1000000)) / (ploidy*avg_cov) # tetraploid size: 3374.929 Mb, haploid size: 843.7321
                message('ploidy for plot is ',stoped_ploidy_ind,'\n')
              }
              # Construct the legend text for portion sizes
              portion_legend <- c(
                paste("Portion-copy ", 1:stoped_ploidy_ind, ": ", round(portion_size[1:stoped_ploidy_ind], digits = 2), " Mb", sep = ""),
                paste("Portion-copy > ", stoped_ploidy_ind, ": ", round(portion_size[stoped_ploidy_ind + 1], digits = 2), " Mb", sep = "")
              )
              legend("topright",
                     pch    = rep(15, 3),
                     col    = c("gray", "gray", "gray","cyan", brewer_palette[2:stoped_ploidy_ind],"gray", alpha("orangered", 0.8)),
                     legend = c(paste("Raw", sep=""),
                                paste("Haploid-peak: ", avg_cov, sep=""),
                                paste("Haploid-fitting:", sep=""),
                                portion_legend,
                                #paste("Portion-copy",1:(stoped_ploidy_ind+1),": ",paste(round(portion_size, digits = 2), sep=","), " Mb", sep=""),
                                paste("Haploid GSE (avg.): ", round(genome_size, digits = 2), " Mb", sep="")
                     ),
                     cex    = 0.8,
                     horiz  = FALSE,
                     box.col="NA")

              dev.off()

              next
            }
            else{
              #abline(v = avg_cov*ploidy_ind, lty=3, col="gray")
              abline(v = meanfit_old*ploidy_ind, lty=3, col="gray")
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
            message('Right now ploidy equal to ploidy_ind\n')
            #lines(histo_raw$V1,  histo_raw$V2, col="blue", lwd=3, lty=2)
            rep_size <- genome_size - sum(portion_size)
            portion_size <- c(portion_size,rep_size)
            if(!ploidy_assign){
              genome_size     = sum(histo_fit$V1 * (histo_fit$V2 / 1000000)) / (ploidy*avg_cov) # tetraploid size: 3374.929 Mb, haploid size: 843.7321
              message('ploidy for plot is ',final_ploidy,'\n')
            }
            # Construct the legend text for portion sizes
            portion_legend <- c(
              paste("Portion-copy ", 1:ploidy, ": ", round(portion_size[1:ploidy], digits = 2), " Mb", sep = ""),
              paste("Portion-copy > ", ploidy, ": ", round(portion_size[ploidy + 1], digits = 2), " Mb", sep = "")
            )
            legend("topright",
                   pch    = rep(15, 3),
                   col    = c("gray", "gray", "gray","cyan", brewer_palette[2:ploidy_ind],"gray", alpha("orangered", 0.8)),
                   legend = c(paste("Raw", sep=""),
                              paste("Haploid-peak: ", avg_cov, sep=""),
                              paste("Haploid-fitting:", sep=""),
                              portion_legend,
                              #paste("Portion-copy",1:(ploidy+1),": ",paste(round(portion_size, digits = 2), sep=","), " Mb", sep=""),
                              paste("Haploid GSE (avg.): ", round(genome_size, digits = 2), " Mb", sep="")
                   ),
                   cex    = 0.8,
                   horiz  = FALSE,
                   box.col="NA")
            #message('ooopppsss.\n')
            dev.off()
          }
        }





        end_time <- initialize_start_time()
        execution_time <- end_time - start_time
        message(paste("Program running time:", sum(execution_time)))
        # grDevices::png file

        grDevices::png(paste(path, output_dir ,sample,"_hap_genome_size_est.png", sep=""),  width = 1200, height = 800, res = 200)
        #grDevices::png(paste(path, output_dir ,sample,"_hap_genome_size_est.grDevices::png", sep=""),  height=4, width=7.08661)
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
            #message("hap peak at", happeak_index, '\n')
            #message("left peak at", range_left, '\n')
            # range_avg_cov <- (happeak_index-range_left):(happeak_index+range_right)
            #avg_cov       <- round(sum(histo_raw[range_avg_cov, 1] * histo_raw[range_avg_cov, 2] / sum(histo_raw[range_avg_cov, 2])), digits = 2)

            # find the valley on the left of het region
            valley_index  <- which.min(histo_raw[1:happeak_index, 2])
            #message(valley_index, '\n')
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
            #abline(v = avg_cov*ploidy_ind, lty=3, col="gray")
            abline(v = meanfit_old*ploidy_ind, lty=3, col="gray")
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
                              portion_legend,
                              #paste("Portion-copy",1:(ploidy_ind+1),": ",paste(round(portion_size, digits = 2), sep=","), " Mb", sep=""),
                              paste("Haploid GSE (avg.): ", round(genome_size, digits = 2), " Mb", sep="")
                   ),
                   cex    = 0.6,
                   horiz  = FALSE,
                   box.col="NA")

            dev.off()
          }
        }
      }
      else{
        message("het_pos is ",het_pos," enter het_pos > 80 process, rescaling... \n")
        if (het_pos >= 80) {
          scaled_factor <- round(het_pos / 60, digits=1)
        }
        scale_flag = TRUE
        #scaled_factor <- round(het_pos / 50)
        message("scaled_factor is ",scaled_factor,"\n")
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
        write.table(file=paste0(path,sample,"_scaled"), x=result_scaled, quote=FALSE, row.names=FALSE, col.names=FALSE)
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
            message('ploidy index : ',ploidy_ind)
            histo_raw <- read.table(paste(path, sample, sep="") ) # from jellyfish
            for(left_fit_ratio in left_fit_ratio_list){
              fit_para_save_list = findGSE_sp(histo = paste(path, sample, sep = ""), sizek = sizek , outdir = paste(path, output_dir, sep = ""), exp_hom = exp_hom_temp,ploidy_ind=ploidy_ind, avg_cov=avg_cov,left_fit_ratio=left_fit_ratio, meanfit_old=meanfit_old, sdfit_old = sdfit_old, scale_flag = scale_flag)
              #grDevices::png(paste("apple_hap_genome_size_est.jpg", sep=""),  height=400, width=708)
              ####
              het_xfit_right_save <- fit_para_save_list[[1]]
              hom_peak_pos_save <- fit_para_save_list[[2]]
              het_peak_pos_save <- fit_para_save_list[[3]]
              meanfit_old= fit_para_save_list[[4]]
              sdfit_old = fit_para_save_list[[5]]
              histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R
              message("rescaling procedure start...\n")
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

              write.table(histo_het_rescale, file = paste0(path, output_dir,sample,"_",ploidy_ind,"_rescale.histo"), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)



              happeak_index_rescale <- which.max(histo_het_rescale$V2) # coverage at het-peak 14, x, ...
              happeak_index_rescale <- histo_het_rescale[happeak_index_rescale,1]
              range_avg_cov_rescale <- (happeak_index_rescale-range_left):(happeak_index_rescale+range_right)
              avg_cov_rescale       <- round(sum(histo_raw_rescale_raw[range_avg_cov_rescale, 1] * histo_raw_rescale_raw[range_avg_cov_rescale, 2] / sum(histo_raw_rescale_raw[range_avg_cov_rescale, 2])), digits = 2)
              message("rescaling procedure end...\n")

              final_fit_peak_pos_ind <- which.max(histo_het_rescale$V2)
              final_fit_peak_pos <-  histo_het_rescale[final_fit_peak_pos_ind,1]
              distance_fit_two_line <- abs(avg_cov_rescale*ploidy_ind - final_fit_peak_pos)
              message('difference of cyan peaks: ',final_fit_peak_pos,', gray peak: ',avg_cov_rescale*ploidy_ind,'\n')
              if(distance_fit_two_line <= 5){
                message(left_fit_ratio,' meet the requirement, continue...\n')
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
            message("hap peak at", happeak_index, '\n')
            #message("left peak at", range_left, '\n')
            if(happeak_index <= range_left){
              range_avg_cov <- 1:(happeak_index+range_right)
            }
            else{
              range_avg_cov <- (happeak_index-range_left):(happeak_index+range_right)
            }

            message('happeak_index is: ',happeak_index,' range_left is ',range_left,'\n')
            avg_cov       <- round(sum(histo_raw[range_avg_cov, 1] * histo_raw[range_avg_cov, 2] / sum(histo_raw[range_avg_cov, 2])), digits = 2)

            # find the valley on the left of het region
            # message('start estimate haploid portion333...: ',happeak_index_rescale,'eee ',length(histo_raw_rescale_raw[1:happeak_index_rescale, 2]),'\n')
            valley_index  <- which.min(histo_raw_rescale_raw[1:happeak_index_rescale, 2])


            #
            # make new histogram with fittings at het-kmers
            histo_fit                    <- histo_raw_rescale_raw

            ## add 0 to histo_fit
            # Find the row number in histo_het_rescale where the first column value is valley_index
            valley_index_k <- which.min(abs(histo_het_rescale[,1] - valley_index))


            message('valley_index is: ',valley_index, ', and valley_index_k is ',valley_index_k,'\n')
            # Calculate the difference
            diff_index <- valley_index - valley_index_k

            # Construct the new vector to assign to histo_fit
            new_values <- c(rep(0, diff_index), histo_het_rescale[1:valley_index_k, 2])
            message('histo_fit the first ',valley_index,' value is: ',new_values,'\n')
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
                 frame.plot = FALSE)
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
              message('max_het_ind is : ',max_het_ind,'\n')
              if (max_het_ind && round(hom_peak_pos_save) < max_het_ind && (max_het_ind - round(hom_peak_pos_save)) > 3 ){
                message('match latter peak pattern, reverse data in ... ploidy ',ploidy_ind,'\n')
                histo_tmp_raw <- histo_tmp
                histo_tmp_raw[target_start_index:length(histo_tmp$V1),2] <- 1
                histo_tmp_raw[target_start_index:target_end_index, 2] <- rev(histo_tmp_raw[round(het_peak_pos_save):round(hom_peak_pos_save), 2])
                write.table(histo_tmp_raw, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
                #lines(histo_tmp_raw$V1,  histo_tmp_raw$V2, col="black", lwd=3)
                continue_reverse_flag = TRUE
              }
              else{
                if(continue_reverse_flag){
                  histo_tmp_raw <- histo_tmp
                  # histo_tmp_raw[which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))):length(histo_tmp_raw$V1),2] <- 1
                  histo_tmp_raw[target_start_index:length(histo_tmp$V1),2] <- 1
                  histo_tmp_raw[target_start_index:target_end_index, 2] <- rev(histo_tmp_raw[round(het_peak_pos_save):round(hom_peak_pos_save), 2])

                  #message('pppploidy_ind range: ',which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))),'\n')
                  write.table(histo_tmp_raw, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)

                }
                else{
                  histo_tmp_raw <- histo_tmp
                  # histo_tmp_raw[which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))):length(histo_tmp_raw$V1),2] <- 1
                  #histo_tmp_raw[target_start_index:length(histo_tmp$V1),2] <- 1
                  #histo_tmp_raw[target_start_index:target_end_index, 2] <- rev(histo_tmp_raw[round(het_peak_pos_save):round(hom_peak_pos_save), 2])

                  #message('pppploidy_ind range: ',which(histo_tmp$V1 == round(avg_cov * (ploidy_ind+0.7))),'\n')
                  write.table(histo_tmp_raw, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
                  #lines(histo_tmp_raw$V1,  histo_tmp_raw$V2, col="black", lwd=3)
                }

              }

            }
            else{
              write.table(histo_tmp, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
            }
            ## end add reverse data part

            # write.table(histo_tmp, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = FALSE)

            exp_hom_temp = avg_cov * (ploidy_ind+1) + avg_cov * 0.1
            histo_raw <- histo_tmp
            message('exp_hom_temp: ',exp_hom_temp,' \n')
            for(left_fit_ratio in left_fit_ratio_list){
              message('Now is fitting left fit ratio: ',left_fit_ratio,'\n')
              fit_para_save_list = findGSE_sp(histo = paste(path, sample,"_",ploidy_ind, sep = ""), sizek = sizek , outdir = paste(path, output_dir, sep = ""), exp_hom = exp_hom_temp,ploidy_ind=ploidy_ind,avg_cov=avg_cov,left_fit_ratio=left_fit_ratio, meanfit_old=meanfit_old, sdfit_old = sdfit_old, scale_flag = scale_flag)
              het_xfit_right_save <- fit_para_save_list[[1]]
              hom_peak_pos_save <- fit_para_save_list[[2]]
              het_peak_pos_save <- fit_para_save_list[[3]]
              histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,"_",ploidy_ind,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R
              # het fitting rescale
              #lines(histo_het$V1,  histo_het$V2, col="deepskyblue", lwd=3, lty=2)
              message('Now scale factor is : ',scaled_factor,'\n')
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
              write.table(histo_het_rescale, file = paste0(path, output_dir,sample,"_",ploidy_ind,"_rescale.histo"), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
              all_slope_data <- abs(diff(diff(histo_het$V2))[(round(avg_cov)*ploidy_ind-10):(round(avg_cov)*ploidy_ind+10)])
              if(!all(all_slope_data > 100)){
                #lines(histo_het_rescale$V1,  histo_het_rescale$V2, col="deepskyblue", lwd=3, lty=2)

                message('index is :',ploidy_ind,' slope is ',all_slope_data)
                success_file <- FALSE
                stoped_ploidy_ind = ploidy_ind - 1
                break
              }
              final_fit_peak_pos_ind <- which.max(histo_het_rescale$V2)

              final_fit_peak_pos <-  histo_het_rescale[final_fit_peak_pos_ind,1]
              distance_fit_two_line <- abs(avg_cov_rescale*ploidy_ind - final_fit_peak_pos)
              message('difference of cyan peaks: ',final_fit_peak_pos,', gray peak: ',avg_cov_rescale*ploidy_ind,'\n')
              if(distance_fit_two_line <= 5){
                message(left_fit_ratio,' meet the requirement, continue...\n')
                break
              }
              message('Do Not meet requirement, try another parameter...\n')
            }

            if(!success_file){

              message('Right now ploidy equal to ploidy_ind\n')

              abline(v = avg_cov_rescale*stoped_ploidy_ind, lty=3, col="gray")

              rep_size <- genome_size - sum(portion_size)
              portion_size <- c(portion_size,rep_size)
              portion_legend <- c(
                paste("Portion-copy ", 1:stoped_ploidy_ind, ": ", round(portion_size[1:stoped_ploidy_ind], digits = 2), " Mb", sep = ""),
                paste("Portion-copy > ", stoped_ploidy_ind, ": ", round(portion_size[stoped_ploidy_ind + 1], digits = 2), " Mb", sep = "")
              )
              legend("topright",
                     pch    = rep(15, 3),
                     col    = c("gray", "gray", "gray","cyan", brewer_palette[2:stoped_ploidy_ind],"gray", alpha("orangered", 0.8)),
                     legend = c(paste("Raw", sep=""),
                                paste("Haploid-peak: ", avg_cov_rescale, sep=""),
                                paste("Haploid-fitting:", sep=""),
                                portion_legend,
                                #paste("Portion-copy",1:(stoped_ploidy_ind+1),": ",paste(round(portion_size, digits = 2), sep=","), " Mb", sep=""),
                                paste("Haploid GSE (avg.): ", round(genome_size, digits = 2), " Mb", sep="")
                     ),
                     cex    = 0.8,
                     horiz  = FALSE,
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
            portion_legend <- c(
              paste("Portion-copy ", 1:ploidy, ": ", round(portion_size[1:ploidy], digits = 2), " Mb", sep = ""),
              paste("Portion-copy > ", ploidy, ": ", round(portion_size[ploidy + 1], digits = 2), " Mb", sep = "")
            )
            legend("topright",
                   pch    = rep(15, 3),
                   col    = c("gray", "gray", "gray","cyan", brewer_palette[2:ploidy_ind],"gray", alpha("orangered", 0.8)),
                   legend = c(paste("Raw", sep=""),
                              paste("Haploid-peak: ", avg_cov_rescale, sep=""),
                              paste("Haploid-fitting", sep=""),
                              portion_legend,
                              #paste("Portion-copy",1:(ploidy_ind+1),": ",paste(round(portion_size, digits = 2), sep=","), " Mb", sep=""),
                              paste("Haploid GSE (avg.): ", round(genome_size, digits = 2), " Mb", sep="")
                   ),
                   cex    = 0.8,
                   horiz  = FALSE,
                   box.col="NA")
            dev.off()
          }
        }

        ## end main code

        end_time <- initialize_start_time()
        execution_time <- end_time - start_time
        message(paste("Program running time:", sum(execution_time)))
        ## grDevices::png file draw


        grDevices::png(paste(path, output_dir ,sub("_scaled$", "", sample),"_hap_genome_size_est.png", sep=""),  width = 1200, height = 800, res = 200)
        par(family = "Helvetica")
        if(!success_file){
          ploidy_raw = ploidy
          ploidy = stoped_ploidy_ind
        }
        for (ploidy_ind in seq(ploidy)){
          if(ploidy_ind == 1){

            avg_cov = 0
            message('ploidy index : ',ploidy_ind)
            #grDevices::png(paste("apple_hap_genome_size_est.jpg", sep=""),  height=400, width=708)
            ####
            histo_raw <- read.table(paste(path, sub("_scaled$", "", sample), sep="") ) # from jellyfish
            histo_het <- read.table(paste(path, output_dir, "v1.94.est.",sample,".genome.size.estimated.k",sizek, "to", sizek, ".fitted_hetfit_count.txt", sep="") ) # from prepare_findGSE.R
            ## re-scale back
            message("rescaling procedure start...\n")
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
            message("rescaling procedure end...\n")

            ## end re-scale back




            # find the haplotype average kmer coverage.
            happeak_index <- which.max(histo_het$V2) # coverage at het-peak 14, x, ...
            #message("hap peak at", happeak_index, '\n')
            #message("left peak at", range_left, '\n')
            range_avg_cov <- (happeak_index-range_left):(happeak_index+range_right)
            avg_cov       <- round(sum(histo_raw[range_avg_cov, 1] * histo_raw[range_avg_cov, 2] / sum(histo_raw[range_avg_cov, 2])), digits = 2)

            # find the valley on the left of het region
            valley_index  <- which.min(histo_raw_rescale_raw[1:happeak_index_rescale, 2])

            #message(valley_index, '\n')
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
                              portion_legend,
                              paste("Portion-copy",1:(ploidy_ind+1),": ",paste(round(portion_size, digits = 2), sep=","), " Mb", sep=""),
                              paste("Haploid GSE (avg.): ", round(genome_size, digits = 2), " Mb", sep="")
                   ),
                   cex    = 0.6,
                   horiz  = FALSE,
                   box.col="NA")
            dev.off()
          }
        }
      }
    }






    #    # grDevices::png file

    #   grDevices::png(paste(path, output_dir ,sample,"_hap_genome_size_est.grDevices::png", sep=""),  width = 800, height = 600)

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
    #       #message("hap peak at", happeak_index, '\n')
    #       #message("left peak at", range_left, '\n')
    #       range_avg_cov <- (happeak_index-range_left):(happeak_index+range_right)
    #       avg_cov       <- round(sum(histo_raw[range_avg_cov, 1] * histo_raw[range_avg_cov, 2] / sum(histo_raw[range_avg_cov, 2])), digits = 2)
    #       # find the valley on the left of het region
    #       valley_index  <- which.min(histo_raw[1:happeak_index, 2])
    #       #message(valley_index, '\n')
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
    #           frame.plot = FALSE)
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
    #       write.table(histo_tmp, file = paste0(path, sample,"_",ploidy_ind), sep = ' ', col.names = F, row.names = F, quote = FALSE)

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
    write.csv(haploid_size_data, file = paste(path, output_dir ,sample,"_haploid_size.csv", sep=""),row.names=FALSE, quote=FALSE)

  }
}
