## Author: Maile Neel
## 2024

## Usage --

## This script is set up to use within the  RProject Potomac.Manuscript.2024.R.Project.Rproj. It tests IBD using Mantel Test of correlations between distances among sites that were calculated in 09A_Calculate_Site_Water_Distances.r and genetic distance calculated using one of the genpop object created in the script 01_Create_adegenet_geninds.R that was written to ./processed.date/list.of.genpop.objects.rds


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Install and Load Packages  ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require(pacman)) install.packages("pacman")
pacman::p_load(adegenet,
               ade4,
               tidyverse,
               here,
               magrittr,
               update = FALSE)


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Read in Required Files  ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## The code below requires the gen.pop objects created in the adegenet 2018_PR_Manuscript_adegenet_create_dataframes.R script.  If the object below is not in the processed.data folder, run script first to generate the object.

list.of.genpop.objects <- readRDS(here("processed.data","list.of.genpop.objects.rds"))

## It also requires the matrices with distances among sites created using the scripts 09a_Calculated_Site_Water_Distances.R.  If the rds files containing these matrices are not in the processed.data folder, that script needs to be run first.

Site.Final.dist <- readRDS(here("processed.data", "Full_PR_Site_Cost_and_Euclidean_Distances_combined_Matrix.rds"))

Site.Final.dist.tidal <- readRDS(here("processed.data", "Tidal_PR_Site_Cost_and_Euclidean_Distances_combined_Matrix.rds"))

Site.Final.dist.nontidal <- readRDS(here("processed.data", "Nontidal_PR_Site_Cost_and_Euclidean_Distances_combined_Matrix.rds"))

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## IBD Tests For All Samples ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

### All Sites ----
#### Calculate genetic distance ----

## The genetic distance options in adegenet are as follows:
## 1) Nei's; 2) Edward's chord distance; #3) Coancestralilty; 4) classic Euclidean; 5) Absolute (Provesti's)

## Using the gen.pop objects in the list read in above, I used method 2 -Angular distance or Edwards' Chord distance - because it is a Euclidean distance measure.

Distgen_all <- adegenet::dist.genpop(list.of.genpop.objects$genpop.allpotomac,
                           method = 2, upper = TRUE)

#### Test IBD ----
 
## Using the distance matrix from the gdistance cost analysis
## Need to make sure the output matrix from gdistance is a dist object (using as.dist)

ibd_allsamp_allsites_cost <- ade4::mantel.randtest(Distgen_all,
                                                   as.dist(Site.Final.dist))

ibd_allsamp_allsites_cost
plot(ibd_allsamp_allsites_cost)

####Basic IBD plot ----
plot(as.dist(Site.Final.dist), Distgen_all, 
     xlab = "Cost Distance (m)", 
     ylab = "Edward's Chord Genetic Distance")
abline(lm(as.vector(Distgen_all) ~as.vector(Site.Final.dist)), 
       col = "black",lty = 1, lwd = 2)

### NonTidal Sites ----

#### Calculate genetic distance ----

Distgen_all_NonTidal <- adegenet::dist.genpop(list.of.genpop.objects$genpop.allpotomac.NonTidal,
                                    method=2)

#### Test IBD ----

ibd_allsamp_NonTidal <- ade4::mantel.randtest(Distgen_all_NonTidal,
                                              as.dist(Site.Final.dist.nontidal))

ibd_allsamp_NonTidal
plot(ibd_allsamp_NonTidal)

####Basic IBD plot ----

plot(as.dist(Site.Final.dist.nontidal), Distgen_all_NonTidal, 
     xlab = "Cost Distance (m)", 
     ylab = "Edward's Chord Distance")
abline(lm(as.vector(Distgen_all_NonTidal) ~as.vector(Site.Final.dist.nontidal)), 
       col = "black",lty = 1, lwd = 2)

###Tidal Sites ----

#### Calculate Genetic distance ----
Distgen_all_Tidal <- adegenet::dist.genpop(list.of.genpop.objects$genpop.allpotomac.Tidal,
                                 method = 2)

#### Test IBD ----

ibd_allsamp_Tidal <- ade4::mantel.randtest(Distgen_all_Tidal,
                                           as.dist(Site.Final.dist.tidal))

ibd_allsamp_Tidal
plot(ibd_allsamp_Tidal)

#### Basic IBD plot ----
plot(as.dist(Site.Final.dist.tidal), Distgen_all_Tidal, 
     xlab = "Cost Distance (m)", 
     ylab = "Edward's Chord Distance")
abline(lm(as.vector(Distgen_all_Tidal) ~as.vector(Site.Final.dist.tidal)), 
       col = "black",lty = 1, lwd = 2)


#+++++++++++++++++++++++++++++++++++++++++++++++++
## IBD Tests For Each MLG In Each Population####
#+++++++++++++++++++++++++++++++++++++++++++++++++

###All Sites ----

#### Calculate Genetic distance ----

Distgen_noreps <- adegenet::dist.genpop(list.of.genpop.objects$genpop.noreps.potomac,
                                        method=2)

#### Test IBD ----

#This is using the distance matrix from the gdistance cost analysis in which the distance is the greater of either cost or Euclidean.

ibd_noreps_allsites_cost <- ade4::mantel.randtest(Distgen_noreps,
                                                  as.dist(Site.Final.dist))

ibd_noreps_allsites_cost
plot(ibd_noreps_allsites_cost)

####Basic IBD plot ----

plot(as.dist(Site.Final.dist), Distgen_noreps)
abline(lm(as.vector(Distgen_noreps)~as.vector(Site.Final.dist)), 
       col="red",lty=2)


### NonTidal Sites only ----

#### Calculate Genetic distance ----

Distgen_noreps_NonTidal <- 
  adegenet::dist.genpop(list.of.genpop.objects$genpop.noreps.potomac.NonTidal,
                                     method = 2)

#### Test IBD ----

ibd_noreps_NonTidal <- ade4::mantel.randtest(Distgen_noreps_NonTidal,
                                             as.dist(Site.Final.dist.nontidal))

plot(ibd_noreps_NonTidal)

#### Basic IBD plot ----
plot(as.dist(Site.Final.dist.nontidal), Distgen_noreps_NonTidal)
abline(lm(as.vector(Distgen_noreps_NonTidal) ~as.vector(Site.Final.dist.nontidal)), 
       col = "red",lty = 2)

### Tidal Sites only ----

#### Calculate Genetic distance ----

Distgen_noreps_Tidal <- adegenet::dist.genpop(list.of.genpop.objects$genpop.noreps.potomac.Tidal,
                                  method=2)

#### Test IBD ----

ibd_noreps_Tidal <- ade4::mantel.randtest(Distgen_noreps_Tidal,
                                          as.dist(Site.Final.dist.tidal),nrepet=999)
ibd_noreps_Tidal
plot(ibd_noreps_Tidal)

#### Basic IBD plot ----

plot(as.dist(Site.Final.dist.tidal), Distgen_noreps_Tidal)
abline(lm(as.vector(Distgen_noreps_Tidal) ~as.vector(Site.Final.dist.tidal)), 
       col = "red",lty = 2)

###Gather up all the IBD tests into a list and pull out the 
###observed correlation  and p values.
Mantel.results <- mget(ls(pattern = "ibd_"))

Mantel.results.summary<-data.frame(Mantel.results.summary = 
                        sapply(Mantel.results, function(x) x$obs), 
                        p.value = sapply(Mantel.results, function(x) x$pvalue))

write.csv(Mantel.results.summary, file = "./processed.data/Potomac.Mantel.Test.Summary.csv")



