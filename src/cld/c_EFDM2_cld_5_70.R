# EFDM : diameter classe model
# Dordogne-Garonne domain (study 752)
# Author: Thierry Bélouard, March 2013
#
#    Copyright (C) 2013  European Union
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, version 2 of the License.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


rm(list = ls())
library(abind)
# install.packages('RODBC')  # Only the first time
library(RODBC)

load("EFDM1_cld_5_70.RData")      # Initial function and variables for Dordogne-Garonne domain

# Data base connection
dbaquitaine <- odbcConnect('EForestSC10')

# Number of steps

# nrofsteps <- 2    # Test
nrofsteps <- 17   # Simulation from 2008 until 2025

# Study domain and scenario of silviculture

simulation.query <- "select * from [DE et scénario] where etude = 752 and plan_libelle = 'gp_prod-mort'"
# Dordogne-Garonne domain and no harvest scenario i.e. production + mortality only
simulation <- sqlQuery(channel = dbaquitaine, query = simulation.query)
nb.simulation <- dim(simulation)[1]

for(num.simul in 1:nb.simulation) {

# Reading of the transition matrix (growth, mortality and harvest all together)

transition.matrix.query <- paste0("select mode_modele, Proba, TxPrelev from [Matrice transition1] where etude = ", 752, " and de = ", simulation$de[num.simul], " and plan_libelle = 'gp_prod-mort' order by mode_modele")
transition.matrix0 <- sqlQuery(channel = dbaquitaine, query = transition.matrix.query)
NumCLD <- dim(transition.matrix0)[1]
GMH <- diag((1-transition.matrix0$TxPrelev)*c(1-transition.matrix0$Proba[-NumCLD], 1))
for(i in 1:(NumCLD-1)) GMH[i+1, i] <- transition.matrix0$Proba[i]*(1-transition.matrix0$TxPrelev[i])
transmats$GMH <- GMH


# Reading of the recruitement

recruit.query <- paste0("select TigesRecrutées from Recrutement1 where ETUDE = 752 and DE = ", simulation$de[num.simul], " and PLAN_LIB = 'gp_prod-mort'")
Recruit <- sqlQuery(channel = dbaquitaine, query = recruit.query)$TigesRecrutées[1]

# Reading of the initial state

#setting up the initial state and the probabilities of activities

initial.state.query <- paste0("select mode_modele as cld, ntig as StemNumber from [état initial1] where etude = 752 and de = ", simulation$de[num.simul], " and plan_libelle = 'gp_prod-mort' order by mode_modele")
statespace <- sqlQuery(channel = dbaquitaine, query = initial.state.query)

statespace <- cbind(statespace, rep(1, dim(statespace)[1]))
dimnames(statespace)[[2]][dim(statespace)[2]] <- "GMH"
for(fname in factnames)
  statespace[[fname]]<-factor(statespace[[fname]],
                      levels = factlvls[[fname]],ordered=TRUE)
statespace <- statespace[do.call(order,statespace[rev(factnames)]),]

initstate <- statespace[[setdiff(names(statespace),union(factnames,actnames))]]
dim(initstate) <- factdims
dimnames(initstate) <- factlvls

actproblist<-list()
for(i in 1:nact) {
  actproblist[[i]]<-array(statespace[[actnames[i]]],dim=factdims)
  dimnames(actproblist[[i]])<-factlvls
}

#finally, the actual simulation

resultstates <- statespace[factnames]
resultstates[["step0"]] <- as.vector(initstate)

state <- initstate
for (i in 1:nrofsteps) {
  state <- newstate(state)
  state[[1]] <- state[[1]] + Recruit
  resultstates[[paste0("step",i)]] <- as.vector(state)
}


#if all went well, the object resultstates now has the state of the population
#at the initial step and further timesteps asked


delete.result.query <- paste0("delete * from EFDMRessource where Etude = 752 and DE = ", simulation[num.simul, ]$de, " and Plan = ", simulation[num.simul, ]$plan)
sqlQuery(channel = dbaquitaine, query = delete.result.query)

for (i in 1:nrofsteps) for (j in 1:13) {
insert.result.query <- paste0("insert into EFDMRessource values (752, ", simulation[num.simul, ]$de, ", ", simulation[num.simul, ]$plan, ", ", 2008 + i, ", ", resultstates[j, ]$cld, ", ", resultstates[j, i + 2], ")")
sqlQuery(channel = dbaquitaine, query = insert.result.query)
}


} # End of the loop : please wait a little

odbcClose(dbaquitaine)