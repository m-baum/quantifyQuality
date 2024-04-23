# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Replication script for                                                      #
# Quantifying the Quality of Configurational Causal Models                    #
#                                                                             #
# Baumgartner & Falk                                                          #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


library(cna)
library(stringr)
library(stringi)
library(useful)
library(QCApro)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("funcs_quality.R")


# \Delta_1, models (4) and (5)
# ----------------------------
### LCR
d1 <- "A*b + c*D <-> E"
m4 <- "A*B + D <-> E"
m5 <- "A*B*D <-> E"

LCR(m4, d1)
LCR(m5, d1)
# ==> Both m4 and m5 are incorrect according to LCR.
# Compare to:
m41 <- "A*b + D <-> E"
m51 <- "A + D <-> E"
LCR(m41, d1)
LCR(m51, d1)
# ==> Both m41 and m51 are correct according to LCR.

### NCR
# Table 1
# Generate submodels of m4
sub.m4 <- get_submodels(m4)
sub.m4
# Generate submodels of m5 
sub.m5 <- get_submodels(m5)
sub.m5
# submodel of d1?
subcheck.m4 <- is.submodel(sub.m4,d1)
subcheck.m5 <- is.submodel(sub.m5,d1)

Table1 <- data.frame(m4=sub.m4, sub_m4_d1=subcheck.m4, m5=sub.m5, sub_m5_d1=subcheck.m5,
                     row.names = NULL)
Table1
# NCR scores
NCR(m4,d1)
NCR(m5,d1)

# \Delta_2, models (4) and (5)
# ----------------------------
d2 <- "A*b*D*F + a*B*C*D <-> E"
m6 <- "A + C + D <-> E"
m7 <- "A*b + C + D <-> E"
m8 <- "A*b*F + C + D <-> E" 

# NCR scores of m6, m7, m8
NCR(m6,d2)
NCR(m7,d2)
NCR(m8,d2)

# Generate submodels of m6, m7, m8
sub.m6 <- get_submodels(m6)
sub.m6
sub.m7 <- get_submodels(m7)
sub.m7
sub.m8 <- get_submodels(m8)
sub.m8

# submodel of d2?
subcheck.m6 <- is.submodel(sub.m6,d2)
subcheck.m7 <- is.submodel(sub.m7,d2)
subcheck.m8 <- is.submodel(sub.m8,d2)

# Table 2
Table2 <- data.frame(m7=sub.m7,sub.m7_sub.m6=is.submodel(sub.m7,m6),
           sub.d2=subcheck.m7, row.names = NULL)
Table2 <- Table2[order(Table2$sub.m7_sub.m6, decreasing = T),]
Table2 <- Table2[order(Table2$sub.d2, decreasing = T),]
row.names(Table2) <- 1:nrow(Table2)
Table2

# Causal expositions of \Delta_3 and (16)
# ---------------------------------------
d3 <- "(A+B*F<->D)*(C+B*f<->E)*(D+E<->G)"
relevancies(d3) # exposition of d3
m16 <- "(A*B<->D)*(D + B*C<->G)"
relevancies(m16)

# Quantifying model qualities, i.e. correctness and completeness
# --------------------------------------------------------------
# Quality of (16) for d3
qual16 <- quality(m16,d3)
qual16
# Breakdown correctness and completeness
qual16[[1]]$correctness
qual16[[1]]$completeness

# Quality of models m1 to m9 for d3 in Tables 6 and 8
T8m1 <- "(A*B <->D)*(B*C <-> G)"
T8m2 <- "(A + B <-> D)*(B*C <-> G)"
T8m3 <- "(A + B*F <-> D)*(B + C <-> G)"
T8m4 <- "A + B <-> G"
T8m5 <- "A + B + E <-> G"
T8m6 <- "A + B + E + D <-> G"
T8m7 <- "A + B + E + D + F <-> G"
T8m8 <- "(A+ B*F <-> D)*(B*f + C <-> E)*(D + E <-> G)"
T8m9 <- "(A+ B*F + H <-> D)*(B*f + C + K <-> E)*(D + E <-> G)"
qualTab8 <- quality(c(T8m1,T8m2,T8m3,T8m4,T8m5,T8m6,T8m7,T8m8,T8m9), d3)
qualTab8
# Breakdown correctness and completeness for model m3
qualTab8[[3]]$correctness
qualTab8[[3]]$completeness
# ==> Change the list index accordingly to get breakdowns for
# the other models.

# Quality of models m6 to m8 for d2 in Table 6b
qualTab6b <- quality(c(m6,m7,m8), d2)
qualTab6b

# Aggregating correctness and completeness
# Qualities in Table 9
# beta = 0.5
qualTab9_beta05 <- quality(c(m16,T8m1,T8m2,T8m3,T8m4,T8m5,T8m6,T8m7,T8m8,T8m9), d3, beta = 0.5)
qualTab9_beta05
qualTab9_beta2 <- quality(c(m16,T8m1,T8m2,T8m3,T8m4,T8m5,T8m6,T8m7,T8m8,T8m9), d3, beta = 2)
qualTab9_beta2
