# calculate the OR for the vaccinee-only test-negative design under
# the null hypothesis of no waning
# OR calculated at the end should be 1 under the null of no waning,
# as the late and early vaccinees should both have (1-VE) the flu
# incidence of unvaccinated persons, and assuming flu vaccination
# doesn't affect nonflu incidence, nonflu should provide a constant
# denominator of incidence.

# Important assumptions as written:
# 1. nonflu unaffected by the vaccine
# 2. nonflu not immunizing -- people remain at risk
# 3. all compared groups are vaccinated after flu season starts,
#    so there is loss from the at-risk group due to flu testing
#    before vaccination, which is differential for the two groups.
# 4. high-risk group for flu incidence is NOT high-risk for nonflu
#    incidence (can change this, and this shifts the bias)

VE <- 0.6
RR <- 1-VE #flu incidence rate ratio for vacc group

# Define number at baseline in high and low risk groups for
# each vaccination treatment (early and late);
# must be equal for early and late by no unmeas confounding
Nb<-c(0,0)
Nb[1]<-100000
Nb[2]<-20000

r<-c(0,0)
r[1]<-.001
r[2]<- r[1]
# nonflu incidence per day in two groups (low and high flu incidence).
# May be unequal among the groups and this changes results

# Specify flu incidence per day in two risk groups
f<-c(0,0)
f[1]<-.001
f[2]<-.004

# day of vaccination for early and late groups, relative to flu season start
d<-c(0,0)
d[1]<-0
d[2]<-60

# day on which incidence is compared between early and late vaccinees
dt<-75

elig=matrix(c(0,0,0,0),nrow=2,ncol=2)
a=elig
incflu = elig
incnon=elig

for (i in (1:2)){    #risk group 1=low, 2=hi
  for (j in (1:2))   #vaccination time 1=early, 2=late
  {elig[i,j] <-Nb[i]*exp(-d[j]*(f[i]+r[i]))        # still eligible for inclusion in study at time of vaccination because no prior flu test
  a[i,j]<-elig[i,j]*exp(-(dt-d[j])*(f[i]*RR+r[i])) # at risk on day dt because no flu test between vax and day t
  incflu[i,j]<-a[i,j]*f[i]*RR                      # incident flu cases on day dt
  incnon[i,j]<-a[i,j]*r[i]                         # incident nonflu cases on day dt
  }}

incnontot<-c(0,0)
incflutot<-c(0,0)
for (j in 1:2) {incflutot[j]<-incflu[1,j]+incflu[2,j]
incnontot[j]<-incnon[1,j]+incnon[2,j]}
OR<-incflutot[1]*incnontot[2]/incflutot[2]/incnontot[1] #comparing early vs. late vaccination; early worse = waning = OR>1
OR

