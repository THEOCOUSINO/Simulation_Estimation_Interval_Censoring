require(VAM)

# paramere Weibull
alpha1<-8e-07
beta1<-2

# parametre Age virtuel
rhop<-0.25
rhopc<-0.75

alpha_ini<-1e-06
beta_ini<-2.5
rho_p_ini<-0.5
rho_pc_ini<-0.5


Cens_h<-(10)*365*24 # date de censure en heure 
Cens_j<-(10)*365 # date de censure en jour
m<-ceiling(Cens_h/8760) # ~ nombre de PM
n<-50 # nb de systemes par jeu de donnees

nb_monte_carlo<-100

(1/(alpha1))^(1/(beta1))

mod_<-NULL

for(i in 1:n){
  
  PM<-ceiling((8760*(1:m)+runif(m,-2100,2100))/24) # les dates de PM
  # les n systemes d'un meme jeu de donnes n'ont pas exactement les meme dates de PM
  # par contre 2 jeux de donnees simules differents ont les meme dates de PM pour les differents systemes
  
  
  mod_[[i]]<-sim.vam(~(ABAO()|Weibull(alpha1,beta1))&(ARA1(rhop)+ARAInf(rhopc)|AtTimes(PM,differentTypeIfCM=TRUE))) # MP - MC
}

Verification_MTTF<-function(donnees){
  MTTF<-c()
  for (i in 1:n){
    MTTF<-c(MTTF,donnees$Time[(donnees$Type==-1) & (donnees$System==i)][1])
  }
  MTTF<-MTTF[!is.na(MTTF)]
  MTTF<-MTTF/360;
  return(MTTF)
}

Verification_pourcentage_defaillance<-function(donnees){
  pourcentage<-c()
  for (i in 1:n){
    pourcentage<-c(pourcentage,length(donnees$Time[(donnees$Type==-1) & (donnees$System==i)])/
                     (length(donnees$Time[(donnees$Type==2) & (donnees$System==i)])+
                        length(donnees$Time[(donnees$Type==0) & (donnees$System==i)])+
                        length(donnees$Time[(donnees$Type==1) & (donnees$System==i)])))
  };
  return(pourcentage)
}

Verification_pourcentage_non_defaillance<-function(donnees){
  pourcentage<-c()
  for (i in 1:n){
    pourcentage<-c(pourcentage,length(donnees$Time[(donnees$Type==1) & (donnees$System==i)])/
                     (length(donnees$Time[(donnees$Type==2) & (donnees$System==i)])+
                        length(donnees$Time[(donnees$Type==0) & (donnees$System==i)])))
  };
  return(pourcentage)
}

Comptage_defaillance<-function(donnees){
  
  Comptage<-data.frame(zero=0*c(1:n),un=0*c(1:n),deux_ou_plus=0*c(1:n))
  Pourcentage<-data.frame(zero=0*c(1:n),un=0*c(1:n),deux_ou_plus=0*c(1:n),nbr_intervalle=0*c(1:n))
  
  for (i in c(1:n)){
    
    Somme_ligne<-0
    donnees_<-data.frame(System=donnees$System[(donnees$System)==i],
                         Time=donnees$Time[(donnees$System)==i],
                         Type=donnees$Type[(donnees$System)==i])
    
    ttCM<-cumsum(donnees_$Type==-1)
    nCM<-ttCM[donnees_$Type>=0]
    
    nombre_defaillance_intervalle<-nCM-c(0,nCM[1:(length(nCM)-1)]) 
    
    for (j in c(1:(length(nombre_defaillance_intervalle)))){
      if(nombre_defaillance_intervalle[j]>0){
        if(nombre_defaillance_intervalle[j]>1){
          Comptage[i,3]<-Comptage[i,3]+1
          Somme_ligne<-Somme_ligne+1
        }
        if((nombre_defaillance_intervalle[j])==1){
          Comptage[i,2]<-Comptage[i,2]+1
          Somme_ligne<-Somme_ligne+1
        }
      }
      if((nombre_defaillance_intervalle[j])==0){
        Comptage[i,1]<-Comptage[i,1]+1
        Somme_ligne<-Somme_ligne+1
      }
    }
    
    Pourcentage[i,4]<-Somme_ligne
    for (k in c(1:3)){
      Pourcentage[i,k]<-(Comptage[i,k])/(Pourcentage[i,4])
    }
    
  }
  return(Pourcentage)
}

l_tilde<-function(donnees,val_ini,model){
  
  
  beta<-val_ini[1]
  rho_p<-val_ini[2]
  rho_pc<-val_ini[3]
  
  Calcul_Vraisemblance_<-0
  Somme_numerateur_total<-0
  Somme_denominateur_total<-0
  
  for (i in 1:n){
    
    
    donnees_actuelles<-data.frame(System=donnees$System[(donnees$System)==i],
                                  Time=donnees$Time[(donnees$System)==i],
                                  Type=donnees$Type[(donnees$System)==i])
    
    PM_i<-c(0,donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
    PM_i<-PM_i[-length(PM_i)]
    PM_i_decallage<-c(donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
    
    Duree_inter_MP<-PM_i_decallage-PM_i
    
    Age_virtuel<-c()
    Age_virtuel<-c(Age_virtuel,0)
    
    
    ttCM<-cumsum(donnees_actuelles$Type==-1)
    nCM<-ttCM[donnees_actuelles$Type>=0]
    
    nombre_defaillance_intervalle<-nCM-c(0,nCM[1:(length(nCM)-1)]) 
    
    if (model=="ARAInfARAInf"){
      
      choix_model<-(((1-rho_p)*(nombre_defaillance_intervalle==0))+((1-rho_pc)*(nombre_defaillance_intervalle>0)))
      
      
      for (j in 1:length(Duree_inter_MP)){
        Age_virtuel<-c(Age_virtuel,(Age_virtuel[j]+Duree_inter_MP[j])*choix_model[j])
      }
      
    }
    
    if (model!="ARAInfARAInf"){
      
      choix_model<-((1-rho_p)*(nombre_defaillance_intervalle==0))*(Duree_inter_MP)+
        ((1-rho_pc)*(nombre_defaillance_intervalle>0))*(Duree_inter_MP)
      
      for (u in 1:length(Duree_inter_MP)){
        Age_virtuel<-c(Age_virtuel,Age_virtuel[u]+
                         ((model=="ARA1ARAInf")*(-rho_pc*Age_virtuel[u])*(nombre_defaillance_intervalle[u]>0))+
                         ((model=="ARAInfARA1")*(-rho_p*Age_virtuel[u])*(nombre_defaillance_intervalle[u]==0))+
                         choix_model[u])
      }
      
    }
    
    
    
    
    Somme_numerateur<-0
    Somme_denominateur<-0
    
    for (l in 1:(length(Duree_inter_MP))){
      Somme_denominateur<-Somme_denominateur+(1*((((Age_virtuel[l]+Duree_inter_MP[l])^(beta))-
                                                    ((Age_virtuel[l])^(beta)))))
      
      Somme_numerateur<-Somme_numerateur+((nombre_defaillance_intervalle[l]>0)*1)
      
    }
    
    Somme_numerateur_total<-Somme_numerateur_total+Somme_numerateur
    Somme_denominateur_total<-Somme_denominateur_total+Somme_denominateur
  }
  
  alpha<-Somme_numerateur_total/Somme_denominateur_total
  
  for (i in 1:n){
    
    
    donnees_actuelles<-data.frame(System=donnees$System[(donnees$System)==i],
                                  Time=donnees$Time[(donnees$System)==i],
                                  Type=donnees$Type[(donnees$System)==i])
    
    PM_i<-c(0,donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
    PM_i<-PM_i[-length(PM_i)]
    PM_i_decallage<-c(donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
    
    Duree_inter_MP<-PM_i_decallage-PM_i
    
    Age_virtuel<-c()
    Age_virtuel<-c(Age_virtuel,0)
    
    
    ttCM<-cumsum(donnees_actuelles$Type==-1)
    nCM<-ttCM[donnees_actuelles$Type>=0]
    
    nombre_defaillance_intervalle<-nCM-c(0,nCM[1:(length(nCM)-1)]) 
    
    if (model=="ARAInfARAInf"){
      
      choix_model<-(((1-rho_p)*(nombre_defaillance_intervalle==0))+((1-rho_pc)*(nombre_defaillance_intervalle>0)))
      
      
      for (j in 1:length(Duree_inter_MP)){
        Age_virtuel<-c(Age_virtuel,(Age_virtuel[j]+Duree_inter_MP[j])*choix_model[j])
      }
      
    }
    
    if (model!="ARAInfARAInf"){
      
      choix_model<-((1-rho_p)*(nombre_defaillance_intervalle==0))*(Duree_inter_MP)+
        ((1-rho_pc)*(nombre_defaillance_intervalle>0))*(Duree_inter_MP)
      
      for (u in 1:length(Duree_inter_MP)){
        Age_virtuel<-c(Age_virtuel,Age_virtuel[u]+
                         ((model=="ARA1ARAInf")*(-rho_pc*Age_virtuel[u])*(nombre_defaillance_intervalle[u]>0))+
                         ((model=="ARAInfARA1")*(-rho_p*Age_virtuel[u])*(nombre_defaillance_intervalle[u]==0))+
                         choix_model[u])
      }
      
    }
    
    premier_facteur<-0
    
    for (l in 1:(length(Duree_inter_MP))){
      premier_facteur<-premier_facteur+(1*(nombre_defaillance_intervalle[l]>0))
    }
    premier_facteur<-log(alpha)*premier_facteur
    
    n_i_facto<-c()
    
    for (j in 1:length(nombre_defaillance_intervalle)){
      n_i_factoriel<-1
      if(nombre_defaillance_intervalle[j]!=0){
        for (k in (1:nombre_defaillance_intervalle[j])){
          n_i_factoriel<-n_i_factoriel*k
        }
      }
      n_i_facto<-c(n_i_facto,n_i_factoriel)
    }
    
    
    second_facteur<-0
    
    for (l in 1:(length(Duree_inter_MP))){
      second_facteur<-second_facteur+(1*((nombre_defaillance_intervalle[l]>0)*log(((Age_virtuel[l]+Duree_inter_MP[l])^(beta))-
                                                                                    ((Age_virtuel[l])^(beta)))))
    }
    
    troisieme_facteur<-0
    
    for (l in 1:(length(Duree_inter_MP))){
      troisieme_facteur<-troisieme_facteur+(1*(nombre_defaillance_intervalle[l]>0))
    }
    
    
    
    Calcul_Vraisemblance_<-Calcul_Vraisemblance_+premier_facteur+second_facteur-troisieme_facteur
  };
  return(-Calcul_Vraisemblance_)
}

l_alpha<-function(donnees,val_ini,model,other_paramters){
  
  alpha<-val_ini[1]
  beta<-other_paramters[1]
  rho_p<-other_paramters[2]
  rho_pc<-other_paramters[3]
  
  Calcul_Vraisemblance_<-0
  
  for (i in 1:n){
    
    
    donnees_actuelles<-data.frame(System=donnees$System[(donnees$System)==i],
                                  Time=donnees$Time[(donnees$System)==i],
                                  Type=donnees$Type[(donnees$System)==i])
    
    PM_i<-c(0,donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
    PM_i<-PM_i[-length(PM_i)]
    PM_i_decallage<-c(donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
    
    Duree_inter_MP<-PM_i_decallage-PM_i
    
    Age_virtuel<-c()
    Age_virtuel<-c(Age_virtuel,0)
    
    ttCM<-cumsum(donnees_actuelles$Type==-1)
    nCM<-ttCM[donnees_actuelles$Type>=0]
    
    nombre_defaillance_intervalle<-nCM-c(0,nCM[1:(length(nCM)-1)]) 
    
    if (model=="ARAInfARAInf"){
      
      choix_model<-(((1-rho_p)*(nombre_defaillance_intervalle==0))+((1-rho_pc)*(nombre_defaillance_intervalle>0)))
      
      
      for (j in 1:length(Duree_inter_MP)){
        Age_virtuel<-c(Age_virtuel,(Age_virtuel[j]+Duree_inter_MP[j])*choix_model[j])
      }
      
    }
    
    if (model!="ARAInfARAInf"){
      
      choix_model<-((1-rho_p)*(nombre_defaillance_intervalle==0))*(Duree_inter_MP)+
        ((1-rho_pc)*(nombre_defaillance_intervalle>0))*(Duree_inter_MP)
      
      for (u in 1:length(Duree_inter_MP)){
        Age_virtuel<-c(Age_virtuel,Age_virtuel[u]+
                         ((model=="ARA1ARAInf")*(-rho_pc*Age_virtuel[u])*(nombre_defaillance_intervalle[u]>0))+
                         ((model=="ARAInfARA1")*(-rho_p*Age_virtuel[u])*(nombre_defaillance_intervalle[u]==0))+
                         choix_model[u])
      }
      
    }
    
    premier_facteur<-0
    
    for (l in 1:(length(Duree_inter_MP))){
      premier_facteur<-premier_facteur+((nombre_defaillance_intervalle[l]==0)*(alpha)*(((Age_virtuel[l]+Duree_inter_MP[l])^(beta))-
                                                                                         ((Age_virtuel[l])^(beta))))
    }
    
    second_facteur<-0
    
    for (l in 1:(length(Duree_inter_MP))){
      second_facteur<-second_facteur+((nombre_defaillance_intervalle[l]>0)*(log(1-(exp((-alpha)*(((Age_virtuel[l]+Duree_inter_MP[l])^(beta))-
                                                                                                   ((Age_virtuel[l])^(beta))))))))
    }
    
    
    
    Calcul_Vraisemblance_<-Calcul_Vraisemblance_-premier_facteur+second_facteur
  };
  return(-Calcul_Vraisemblance_)
}

l_alpha_chapeau<-function(donnees_,model_,other_paramters_){
  
  beta<-other_paramters_[1]
  rho_p<-other_paramters_[2]
  rho_pc<-other_paramters_[3]
  
  Somme_numerateur_total<-0
  Somme_denominateur_total<-0
  
  for (i in 1:n){
    
    
    donnees_actuelles<-data.frame(System=donnees_$System[(donnees_$System)==i],
                                  Time=donnees_$Time[(donnees_$System)==i],
                                  Type=donnees_$Type[(donnees_$System)==i])
    
    PM_i<-c(0,donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
    PM_i<-PM_i[-length(PM_i)]
    PM_i_decallage<-c(donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
    
    Duree_inter_MP<-PM_i_decallage-PM_i
    
    Age_virtuel<-c()
    Age_virtuel<-c(Age_virtuel,0)
    
    
    ttCM<-cumsum(donnees_actuelles$Type==-1)
    nCM<-ttCM[donnees_actuelles$Type>=0]
    
    nombre_defaillance_intervalle<-nCM-c(0,nCM[1:(length(nCM)-1)]) 
    
    if (model_=="ARAInfARAInf"){
      
      choix_model<-(((1-rho_p)*(nombre_defaillance_intervalle==0))+((1-rho_pc)*(nombre_defaillance_intervalle>0)))
      
      
      for (m in 1:length(Duree_inter_MP)){
        Age_virtuel<-c(Age_virtuel,(Age_virtuel[m]+Duree_inter_MP[m])*choix_model[m])
      }
      
    }
    
    if (model_!="ARAInfARAInf"){
      
      choix_model<-((1-rho_p)*(nombre_defaillance_intervalle==0))*(Duree_inter_MP)+
        ((1-rho_pc)*(nombre_defaillance_intervalle>0))*(Duree_inter_MP)
      
      for (a in 1:length(Duree_inter_MP)){
        Age_virtuel<-c(Age_virtuel,Age_virtuel[a]+
                         ((model_=="ARA1ARAInf")*(-rho_pc*Age_virtuel[a])*(nombre_defaillance_intervalle[a]>0))+
                         ((model_=="ARAInfARA1")*(-rho_p*Age_virtuel[a])*(nombre_defaillance_intervalle[a]==0))+
                         choix_model[a])
      }
      
    }
    
    Somme_numerateur<-0
    Somme_denominateur<-0
    
    for (l in 1:(length(Duree_inter_MP))){
      Somme_denominateur<-Somme_denominateur+(1*((((Age_virtuel[l]+Duree_inter_MP[l])^(beta))-
                                                    ((Age_virtuel[l])^(beta)))))
      
      Somme_numerateur<-Somme_numerateur+((nombre_defaillance_intervalle[l]>0)*1)
      
    }
    
    Somme_numerateur_total<-Somme_numerateur_total+Somme_numerateur
    Somme_denominateur_total<-Somme_denominateur_total+Somme_denominateur
  }
  
  alpha_tilde<-Somme_numerateur_total/Somme_denominateur_total
  
  Resultat<-optim(alpha_tilde,l_alpha,donnees=donnees_,model="ARA1ARAInf",other_paramters=other_paramters_,lower=c(1e-14),upper=c(1e-04),method="Brent"
  )
  
  if(Resultat$convergence!=0){
    print("3")
  }
  
  alpha<-Resultat$par[1];
  return(alpha)
}


l_trois_parametres<-function(donnees,val_ini,model,other_paramters){
  
  alpha<-l_alpha_chapeau(donnees,model,val_ini)
  
  beta<-val_ini[1]
  rho_p<-val_ini[2]
  rho_pc<-val_ini[3]
  
  
  Calcul_Vraisemblance_<-0
  
  for (i in 1:n){
    
    
    donnees_actuelles<-data.frame(System=donnees$System[(donnees$System)==i],
                                  Time=donnees$Time[(donnees$System)==i],
                                  Type=donnees$Type[(donnees$System)==i])
    
    PM_i<-c(0,donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
    PM_i<-PM_i[-length(PM_i)]
    PM_i_decallage<-c(donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
    
    Duree_inter_MP<-PM_i_decallage-PM_i
    
    Age_virtuel<-c()
    Age_virtuel<-c(Age_virtuel,0)
    
    
    ttCM<-cumsum(donnees_actuelles$Type==-1)
    nCM<-ttCM[donnees_actuelles$Type>=0]
    
    nombre_defaillance_intervalle<-nCM-c(0,nCM[1:(length(nCM)-1)]) 
    
    if (model=="ARAInfARAInf"){
      
      choix_model<-(((1-rho_p)*(nombre_defaillance_intervalle==0))+((1-rho_pc)*(nombre_defaillance_intervalle>0)))
      
      
      for (j in 1:length(Duree_inter_MP)){
        Age_virtuel<-c(Age_virtuel,(Age_virtuel[j]+Duree_inter_MP[j])*choix_model[j])
      }
      
    }
    
    if (model!="ARAInfARAInf"){
      
      choix_model<-((1-rho_p)*(nombre_defaillance_intervalle==0))*(Duree_inter_MP)+
        ((1-rho_pc)*(nombre_defaillance_intervalle>0))*(Duree_inter_MP)
      
      for (u in 1:length(Duree_inter_MP)){
        Age_virtuel<-c(Age_virtuel,Age_virtuel[u]+
                         ((model=="ARA1ARAInf")*(-rho_pc*Age_virtuel[u])*(nombre_defaillance_intervalle[u]>0))+
                         ((model=="ARAInfARA1")*(-rho_p*Age_virtuel[u])*(nombre_defaillance_intervalle[u]==0))+
                         choix_model[u])
      }
      
    }
    
    
    premier_facteur<-0
    
    for (l in 1:(length(Duree_inter_MP))){
      premier_facteur<-premier_facteur+((nombre_defaillance_intervalle[l]==0)*((alpha)*(((Age_virtuel[l]+Duree_inter_MP[l])^(beta))-
                                                                                          ((Age_virtuel[l])^(beta)))))
    }
    
    second_facteur<-0
    
    for (l in 1:(length(Duree_inter_MP))){
      second_facteur<-second_facteur+(nombre_defaillance_intervalle[l]>0)*(log(1-(exp((-alpha)*(((Age_virtuel[l]+Duree_inter_MP[l])^(beta))-
                                                                                                  ((Age_virtuel[l])^(beta)))))))
    }
    
    
    
    Calcul_Vraisemblance_<-Calcul_Vraisemblance_-premier_facteur+second_facteur
  };
  return(-Calcul_Vraisemblance_)
}

Vraisemblance_n_syst_MC<-function(donnees,val_ini,model){
  
  
  beta<-val_ini[1]
  rho_p<-val_ini[2]
  rho_pc<-val_ini[3]
  
  Calcul_Vraisemblance_<-0
  Somme_numerateur_total<-0
  Somme_denominateur_total<-0
  
  for (i in 1:n){
    
    
    donnees_actuelles<-data.frame(System=donnees$System[(donnees$System)==i],
                                  Time=donnees$Time[(donnees$System)==i],
                                  Type=donnees$Type[(donnees$System)==i])
    
    PM_i<-c(0,donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
    PM_i<-PM_i[-length(PM_i)]
    PM_i_decallage<-c(donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
    
    Duree_inter_MP<-PM_i_decallage-PM_i
    
    Age_virtuel<-c()
    Age_virtuel<-c(Age_virtuel,0)
    
    
    ttCM<-cumsum(donnees_actuelles$Type==-1)
    nCM<-ttCM[donnees_actuelles$Type>=0]
    
    nombre_defaillance_intervalle<-nCM-c(0,nCM[1:(length(nCM)-1)]) 
    
    if (model=="ARAInfARAInf"){
      
      choix_model<-(((1-rho_p)*(nombre_defaillance_intervalle==0))+((1-rho_pc)*(nombre_defaillance_intervalle>0)))
      
      
      for (j in 1:length(Duree_inter_MP)){
        Age_virtuel<-c(Age_virtuel,(Age_virtuel[j]+Duree_inter_MP[j])*choix_model[j])
      }
      
    }
    
    if (model!="ARAInfARAInf"){
      
      choix_model<-((1-rho_p)*(nombre_defaillance_intervalle==0))*(Duree_inter_MP)+
        ((1-rho_pc)*(nombre_defaillance_intervalle>0))*(Duree_inter_MP)
      
      
      for (u in 1:length(Duree_inter_MP)){
        Age_virtuel<-c(Age_virtuel,Age_virtuel[u]+
                         ((model=="ARA1ARAInf")*(-rho_pc*Age_virtuel[u])*(nombre_defaillance_intervalle[u]>0))+
                         ((model=="ARAInfARA1")*(-rho_p*Age_virtuel[u])*(nombre_defaillance_intervalle[u]==0))+
                         choix_model[u])
      }
      
    }
    
    
    
    Somme_numerateur<-0
    Somme_denominateur<-0
    
    for (l in 1:(length(Duree_inter_MP))){
      Somme_denominateur<-Somme_denominateur+(((Age_virtuel[l]+Duree_inter_MP[l])^(beta))-
                                                ((Age_virtuel[l])^(beta)))
      
      Somme_numerateur<-Somme_numerateur+nombre_defaillance_intervalle[l]
      
    }
    
    Somme_numerateur_total<-Somme_numerateur_total+Somme_numerateur
    Somme_denominateur_total<-Somme_denominateur_total+Somme_denominateur
  }
  
  alpha<-Somme_numerateur_total/Somme_denominateur_total
  
  for (i in 1:n){
    
    
    donnees_actuelles<-data.frame(System=donnees$System[(donnees$System)==i],
                                  Time=donnees$Time[(donnees$System)==i],
                                  Type=donnees$Type[(donnees$System)==i])
    
    PM_i<-c(0,donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
    PM_i<-PM_i[-length(PM_i)]
    PM_i_decallage<-c(donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
    
    Duree_inter_MP<-PM_i_decallage-PM_i
    
    Age_virtuel<-c()
    Age_virtuel<-c(Age_virtuel,0)
    
    
    ttCM<-cumsum(donnees_actuelles$Type==-1)
    nCM<-ttCM[donnees_actuelles$Type>=0]
    
    nombre_defaillance_intervalle<-nCM-c(0,nCM[1:(length(nCM)-1)]) 
    
    if (model=="ARAInfARAInf"){
      
      choix_model<-(((1-rho_p)*(nombre_defaillance_intervalle==0))+((1-rho_pc)*(nombre_defaillance_intervalle>0)))
      
      
      for (j in 1:length(Duree_inter_MP)){
        Age_virtuel<-c(Age_virtuel,(Age_virtuel[j]+Duree_inter_MP[j])*choix_model[j])
      }
      
    }
    
    if (model!="ARAInfARAInf"){
      
      choix_model<-((1-rho_p)*(nombre_defaillance_intervalle==0))*(Duree_inter_MP)+
        ((1-rho_pc)*(nombre_defaillance_intervalle>0))*(Duree_inter_MP)
      
      for (u in 1:length(Duree_inter_MP)){
        Age_virtuel<-c(Age_virtuel,Age_virtuel[u]+
                         ((model=="ARA1ARAInf")*(-rho_pc*Age_virtuel[u])*(nombre_defaillance_intervalle[u]>0))+
                         ((model=="ARAInfARA1")*(-rho_p*Age_virtuel[u])*(nombre_defaillance_intervalle[u]==0))+
                         choix_model[u])
      }
      
    }
    
    premier_facteur<-0
    
    for (l in 1:(length(Duree_inter_MP))){
      premier_facteur<-premier_facteur+nombre_defaillance_intervalle[l]
    }
    premier_facteur<-log(alpha)*premier_facteur
    
    n_i_facto<-c()
    
    for (j in 1:length(nombre_defaillance_intervalle)){
      n_i_factoriel<-1
      if(nombre_defaillance_intervalle[j]!=0){
        for (k in (1:nombre_defaillance_intervalle[j])){
          n_i_factoriel<-n_i_factoriel*k
        }
      }
      n_i_facto<-c(n_i_facto,n_i_factoriel)
    }
    
    
    second_facteur<-0
    
    for (l in 1:(length(Duree_inter_MP))){
      second_facteur<-second_facteur+((nombre_defaillance_intervalle[l])*log(((Age_virtuel[l]+Duree_inter_MP[l])^(beta))-
                                                                               ((Age_virtuel[l])^(beta))))
    }
    
    troisieme_facteur<-0
    
    for (l in 1:(length(Duree_inter_MP))){
      troisieme_facteur<-troisieme_facteur+nombre_defaillance_intervalle[l]
    }
    
    quatrieme_facteur<-0
    for (l in 1:(length(Duree_inter_MP))){
      quatrieme_facteur<-quatrieme_facteur+(log(n_i_facto[l]))
    }
    
    
    Calcul_Vraisemblance_<-Calcul_Vraisemblance_+premier_facteur+second_facteur-troisieme_facteur-quatrieme_facteur
  };
  return(-Calcul_Vraisemblance_)
}


Estimation_Lecture<-function(model_,nb_monte_carlo_,donnees_globales){
  
  Donnees_sorties<-NULL
  Donnees_sorties_passage<-NULL
  
  Resultat_exact<-NULL
  Resultat_global_connu<-NULL
  Resultat_alpha_connu<-NULL
  Resultat_global_inconnu<-NULL
  Resultat_alpha_inconnu<-NULL
  Resultat_global_inconnu_bis<-NULL
  Resultat_alpha_inconnu_bis<-NULL
  ests_cens_1<-NULL
  ests_cens_2<-NULL
  
  indicateur<-0
  
  Ratio_MP_MPC<-c()
  Defaillances<-c()
  Nombre_MP<-c()
  Nombre_MPC<-c()
  
  MTTF<-c()
  pourcentage_defaillance<-c()
  pourcentage_non_defaillance<-c()
  defaillance<-NULL
  
  for (j in 1:nb_monte_carlo_){
    
    valeurs_initiales_1<-c(beta_ini,rho_p_ini,rho_pc_ini)
    
    print(j)
    
    donnees_<-NULL
    
    for(i in 1:n){
      donnees_<-rbind(donnees_,cbind("System"=i,simulate(mod_[[i]],T>(RC=Cens_j),nb.system = 1)))
    }
    
    
    Ratio_MP_MPC<-c(Ratio_MP_MPC,nrow(donnees_[donnees_$Type==1,])/nrow(donnees_[donnees_$Type==2,]))
    Defaillances<-c(Defaillances,nrow(donnees_[donnees_$Type==-1,]))
    Nombre_MP<-c(Nombre_MP,nrow(donnees_[donnees_$Type==1,]))
    Nombre_MPC<-c(Nombre_MPC,nrow(donnees_[donnees_$Type==2,]))
    
    MTTF<-c(MTTF,Verification_MTTF(donnees_))
    pourcentage_defaillance<-c(pourcentage_defaillance,Verification_pourcentage_defaillance(donnees_))
    pourcentage_non_defaillance<-c(pourcentage_non_defaillance,Verification_pourcentage_non_defaillance(donnees_))
    defaillance<-rbind(defaillance,Comptage_defaillance(donnees_))
    
    methEst_exact<-mle.vam(System & Time & Type~(ABAO()|Weibull(alpha_ini,beta_ini))&(ARA1(rho_p_ini)+ARAInf(rho_pc_ini)),data=donnees_)
    exacts<-run(methEst_exact,lower=c(1,0,0),upper=c(5,1,1),method="L-BFGS-B")
    
    if(exacts$convergence!=0){
      print("1")
    }
    
    Resultat_exact<-rbind(Resultat_exact,exacts$par)
    
    Resultat<-optim(valeurs_initiales_1,Vraisemblance_n_syst_MC,donnees=donnees_,model="ARA1ARAInf",lower=c(1,0,0),upper=c(5,1,1),method="L-BFGS-B")
    
    if(Resultat$convergence!=0){
      print("2")
    }
    
    Resultat_global_connu<-rbind(Resultat_global_connu,Resultat$par)
    
    beta<-Resultat$par[1]
    rho_p<-Resultat$par[2]
    rho_pc<-Resultat$par[3]
    
    Somme_numerateur_total<-0
    Somme_denominateur_total<-0
    
    
    for (i in 1:n){
      
      
      donnees_actuelles<-data.frame(System=donnees_$System[(donnees_$System)==i],
                                    Time=donnees_$Time[(donnees_$System)==i],
                                    Type=donnees_$Type[(donnees_$System)==i])
      
      PM_i<-c(0,donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
      PM_i<-PM_i[-length(PM_i)]
      PM_i_decallage<-c(donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
      
      Duree_inter_MP<-PM_i_decallage-PM_i
      
      Age_virtuel<-c()
      Age_virtuel<-c(Age_virtuel,0)
      
      
      ttCM<-cumsum(donnees_actuelles$Type==-1)
      nCM<-ttCM[donnees_actuelles$Type>=0]
      
      nombre_defaillance_intervalle<-nCM-c(0,nCM[1:(length(nCM)-1)]) 
      
      if (model_=="ARAInfARAInf"){
        
        choix_model<-(((1-rho_p)*(nombre_defaillance_intervalle==0))+((1-rho_pc)*(nombre_defaillance_intervalle>0)))
        
        
        for (p in 1:length(Duree_inter_MP)){
          Age_virtuel<-c(Age_virtuel,(Age_virtuel[p]+Duree_inter_MP[p])*choix_model[p])
        }
        
      }
      
      if (model_!="ARAInfARAInf"){
        
        choix_model<-((1-rho_p)*(nombre_defaillance_intervalle==0))*(Duree_inter_MP)+
          ((1-rho_pc)*(nombre_defaillance_intervalle>0))*(Duree_inter_MP)
        
        for (u in 1:length(Duree_inter_MP)){
          Age_virtuel<-c(Age_virtuel,Age_virtuel[u]+
                           ((model_=="ARA1ARAInf")*(-rho_pc*Age_virtuel[u])*(nombre_defaillance_intervalle[u]>0))+
                           ((model_=="ARAInfARA1")*(-rho_p*Age_virtuel[u])*(nombre_defaillance_intervalle[u]==0))+
                           choix_model[u])
        }
        
      }
      
      Somme_numerateur<-0
      Somme_denominateur<-0
      
      for (l in 1:(length(Duree_inter_MP))){
        Somme_denominateur<-Somme_denominateur+(((Age_virtuel[l]+Duree_inter_MP[l])^(beta))-
                                                  ((Age_virtuel[l])^(beta)))
        
        Somme_numerateur<-Somme_numerateur+nombre_defaillance_intervalle[l]
        
      }
      
      
      Somme_numerateur_total<-Somme_numerateur_total+Somme_numerateur
      Somme_denominateur_total<-Somme_denominateur_total+Somme_denominateur
    }
    
    alpha<-Somme_numerateur_total/Somme_denominateur_total
    
    Resultat_alpha_connu<-rbind(Resultat_alpha_connu,alpha)
    
    Resultat<-optim(valeurs_initiales_1,l_trois_parametres,donnees=donnees_,model="ARA1ARAInf",
                    lower=c(1,0,0),upper=c(5,1,1),method="L-BFGS-B")
    
    if(Resultat$convergence!=0){
      print("3 bis")
    }
    
    beta_2<-Resultat$par[1]
    rho_p_2<-Resultat$par[2]
    rho_pc_2<-Resultat$par[3]
    alpha_2<-l_alpha_chapeau(donnees_,"ARA1ARAInf",c(beta_2,rho_p_2,rho_pc_2))
    
    Resultat_alpha_inconnu<-rbind(Resultat_alpha_inconnu,alpha_2)
    Resultat_global_inconnu<-rbind(Resultat_global_inconnu,c(beta_2,rho_p_2,rho_pc_2))
    
    Resultat<-optim(valeurs_initiales_1,l_tilde,donnees=donnees_,model="ARA1ARAInf",lower=c(1,0,0),upper=c(5,1,1),method="L-BFGS-B")
    
    if(Resultat$convergence!=0){
      print("4")
    }
    
    Resultat_global_inconnu_bis<-rbind(Resultat_global_inconnu_bis,Resultat$par)
    
    beta_1<-Resultat$par[1]
    rho_p_1<-Resultat$par[2]
    rho_pc_1<-Resultat$par[3]
    
    Somme_numerateur_total<-0
    Somme_denominateur_total<-0
    
    for (i in 1:n){
      
      
      donnees_actuelles<-data.frame(System=donnees_$System[(donnees_$System)==i],
                                    Time=donnees_$Time[(donnees_$System)==i],
                                    Type=donnees_$Type[(donnees_$System)==i])
      
      PM_i<-c(0,donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
      PM_i<-PM_i[-length(PM_i)]
      PM_i_decallage<-c(donnees_actuelles$Time[donnees_actuelles$Type>=0]*((donnees_actuelles$Type[donnees_actuelles$Type>=0]==1)|donnees_actuelles$Type[donnees_actuelles$Type>=0]==2|donnees_actuelles$Type[donnees_actuelles$Type>=0]==0))
      
      Duree_inter_MP<-PM_i_decallage-PM_i
      
      Age_virtuel<-c()
      Age_virtuel<-c(Age_virtuel,0)
      
      
      ttCM<-cumsum(donnees_actuelles$Type==-1)
      nCM<-ttCM[donnees_actuelles$Type>=0]
      
      nombre_defaillance_intervalle<-nCM-c(0,nCM[1:(length(nCM)-1)]) 
      
      if (model_=="ARAInfARAInf"){
        
        choix_model<-(((1-rho_p_1)*(nombre_defaillance_intervalle==0))+((1-rho_pc_1)*(nombre_defaillance_intervalle>0)))
        
        
        for (m in 1:length(Duree_inter_MP)){
          Age_virtuel<-c(Age_virtuel,(Age_virtuel[m]+Duree_inter_MP[m])*choix_model[m])
        }
        
      }
      
      if (model_!="ARAInfARAInf"){
        
        choix_model<-((1-rho_p_1)*(nombre_defaillance_intervalle==0))*(Duree_inter_MP)+
          ((1-rho_pc_1)*(nombre_defaillance_intervalle>0))*(Duree_inter_MP)
        
        for (a in 1:length(Duree_inter_MP)){
          Age_virtuel<-c(Age_virtuel,Age_virtuel[a]+
                           ((model_=="ARA1ARAInf")*(-rho_pc_1*Age_virtuel[a])*(nombre_defaillance_intervalle[a]>0))+
                           ((model_=="ARAInfARA1")*(-rho_p_1*Age_virtuel[a])*(nombre_defaillance_intervalle[a]==0))+
                           choix_model[a])
        }
        
      }
      
      Somme_numerateur<-0
      Somme_denominateur<-0
      
      for (l in 1:(length(Duree_inter_MP))){
        Somme_denominateur<-Somme_denominateur+(1*((((Age_virtuel[l]+Duree_inter_MP[l])^(beta_1))-
                                                      ((Age_virtuel[l])^(beta_1)))))
        
        Somme_numerateur<-Somme_numerateur+((nombre_defaillance_intervalle[l]>0)*1)
        
      }
      
      Somme_numerateur_total<-Somme_numerateur_total+Somme_numerateur
      Somme_denominateur_total<-Somme_denominateur_total+Somme_denominateur
    }
    
    alpha_1<-Somme_numerateur_total/Somme_denominateur_total
    
    Resultat_alpha_inconnu_bis<-rbind(Resultat_alpha_inconnu_bis,alpha_1)
    
    
    ttCM<-cumsum(donnees_$Type==-1) 
    nCM<-ttCM[donnees_$Type>=0] 
    ttCMenPM<-rep(0,length(ttCM))
    ttCMenPM[donnees_$Type>=0]<-nCM
    ttCM<-ttCM-cummax(ttCMenPM)
    
    donnees_<-data.frame(donnees_,TimeCens1=donnees_$Time,TimeCens2=donnees_$Time)
    
    # on construit les dates correspondant a l'idee 1 modifiee
    # les bornes inferieures successives de intervalles de censure
    
    PMetCens_prec<-c(0,donnees_$Time[donnees_$Type>=0]*((donnees_$Type[donnees_$Type>=0]==1)|donnees_$Type[donnees_$Type>=0]==2))
    
    # les bornes superieures successives de intervalles de censure
    PMetCens_suiv<-donnees_$Time[donnees_$Type>=0] # le nombre de CM total dans les intervalles de censures successifs
    nbCMentre_prec_et_suiv<-nCM-c(0,nCM[1:(length(nCM)-1)]) # pour chaqu'une des CM successive de donnees a quel intervalle de censure elle correspond
    i_PM<-(cumsum(donnees_$Type>=0))[donnees_$Type==-1]+1
    
    donnees_$TimeCens1[donnees_$Type==-1]<-PMetCens_prec[i_PM]+
      (PMetCens_suiv[i_PM]-PMetCens_prec[i_PM])/
      (nbCMentre_prec_et_suiv[i_PM]+1)*ttCM[donnees_$Type==-1]
    
    # on construit les dates correspondant a l'idee 2
    ind<-ttCM[donnees_$Type==-1]==1
    donnees_$TimeCens2[ttCM==1]<-PMetCens_prec[i_PM[ind]]+
      (PMetCens_suiv[i_PM[ind]]-PMetCens_prec[i_PM[ind]])*runif(sum(ttCM==1))
    k<-2
    while(k<=max(ttCM)){
      ind<-which(ttCM==k)
      donnees_$TimeCens2[ind]<-donnees_$TimeCens2[ind-1]+
        (PMetCens_suiv[i_PM[ttCM[donnees_$Type==-1]==k]]-donnees_$TimeCens2[ind-1])*runif(length(ind))
      k<-k+1
    }
    
    methEst_cens_1<-mle.vam(System & TimeCens1 & Type~(ABAO()|Weibull(alpha_ini,beta_ini))&(ARA1(rho_p_ini)+ARAInf(rho_pc_ini)),data=donnees_)
    methEst_cens_2<-mle.vam(System & TimeCens2 & Type~(ABAO()|Weibull(alpha_ini,beta_ini))&(ARA1(rho_p_ini)+ARAInf(rho_pc_ini)),data=donnees_)
    
    cens_1<-run(methEst_cens_1,lower=c(1,0,0),upper=c(5,1,1),method="L-BFGS-B")
    cens_2<-run(methEst_cens_2,lower=c(1,0,0),upper=c(5,1,1),method="L-BFGS-B")
    
    if(cens_1$convergence!=0){
      print("5")
    }
    
    if(cens_2$convergence!=0){
      print("6")
    }
    
    ests_cens_1<-rbind(ests_cens_1,cens_1$par)
    ests_cens_2<-rbind(ests_cens_2,cens_2$par)
    
    Donnees_sorties_passage<-rbind(Donnees_sorties_passage,donnees_)
    
  }
  
  print(mean(Ratio_MP_MPC))
  print(mean(Defaillances))
  print(mean(Nombre_MP))
  print(mean(Nombre_MPC))
  
  boxplot(MTTF)
  boxplot(pourcentage_defaillance)
  boxplot(pourcentage_non_defaillance)
  print(colMeans(defaillance))
  
  print(c(colMeans(Resultat_alpha_connu),colMeans(Resultat_global_connu)))
  print(c(colMeans(Resultat_alpha_inconnu),colMeans(Resultat_global_inconnu)))
  print(c(colMeans(Resultat_alpha_inconnu_bis),colMeans(Resultat_global_inconnu_bis)))
  print(colMeans(Resultat_exact))
  print(colMeans(ests_cens_1))
  print(colMeans(ests_cens_2))
  
  par(mfrow=c(2,2))
  leg<-c("beta","rho p","rho pc")
  param0<-c(beta1,rhop,rhopc)
  
  boxplot(Resultat_exact[,2],ests_cens_1[,2],ests_cens_2[,2],Resultat_global_connu[,1],Resultat_global_inconnu_bis[,1],
          Resultat_global_inconnu[,1],ylab=leg[1],ylim=c(0,5))
  abline(h=param0[1])
  
  boxplot(Resultat_exact[,3],ests_cens_1[,3],ests_cens_2[,3],Resultat_global_connu[,2],Resultat_global_inconnu_bis[,2],
          Resultat_global_inconnu[,2],ylab=leg[2],ylim=c(0,1))
  abline(h=param0[2])
  
  boxplot(Resultat_exact[,4],ests_cens_1[,4],ests_cens_2[,4],Resultat_global_connu[,3],Resultat_global_inconnu_bis[,3],
          Resultat_global_inconnu[,3],ylab=leg[3],ylim=c(0,1))
  abline(h=param0[3])
  
  eta_exacts<-c()
  eta_cens_connu<-c()
  eta_cens_inconnu<-c()
  eta_cens_inconnu_bis<-c()
  eta_ests_cens_1<-c()
  eta_ests_cens_2<-c()
  
  for(j in 1:nb_monte_carlo){
    eta_exacts<-c(eta_exacts,(1/(Resultat_exact[j,1]))^(1/(Resultat_exact[j,2])))
    eta_cens_connu<-c(eta_cens_connu,(1/(Resultat_alpha_connu[j]))^(1/(Resultat_global_connu[j,1])))
    eta_cens_inconnu<-c(eta_cens_inconnu,(1/(Resultat_alpha_inconnu[j]))^(1/(Resultat_global_inconnu[j,1])))
    eta_cens_inconnu_bis<-c(eta_cens_inconnu_bis,(1/(Resultat_alpha_inconnu_bis[j]))^(1/(Resultat_global_inconnu_bis[j,1])))
    eta_ests_cens_1<-c(eta_ests_cens_1,(1/(ests_cens_1[j,1]))^(1/(ests_cens_1[j,2])))
    eta_ests_cens_2<-c(eta_ests_cens_2,(1/(ests_cens_2[j,1]))^(1/(ests_cens_2[j,2])))
  }
  
  
  leg<-c("eta")
  param0<-c((1/(alpha1))^(1/(beta1)))
  boxplot(eta_exacts,eta_ests_cens_1,eta_ests_cens_2,eta_cens_connu,eta_cens_inconnu,eta_cens_inconnu_bis,ylab="eta")
  abline(h=param0[1])
  
  print(mean(eta_cens_connu))
  print(mean(eta_cens_inconnu))
  print(mean(eta_cens_inconnu_bis))
  print(mean(eta_exacts))
  print(mean(eta_ests_cens_1))
  print(mean(eta_ests_cens_2))
  
  Donnees_sorties[[1]]<-Donnees_sorties_passage
  Donnees_sorties[[2]]<-cbind(Resultat_alpha_connu,Resultat_global_connu)
  Donnees_sorties[[3]]<-cbind(Resultat_alpha_inconnu,Resultat_global_inconnu)
  Donnees_sorties[[4]]<-cbind(Resultat_alpha_inconnu_bis,Resultat_global_inconnu_bis)
  Donnees_sorties[[5]]<-Resultat_exact
  Donnees_sorties[[6]]<-ests_cens_1
  Donnees_sorties[[7]]<-ests_cens_2;
  
  return(Donnees_sorties)
  
  
}

Donnees_sorties<-Estimation_Lecture("ARA1ARAInf",nb_monte_carlo,Echantillon_global)


simARAInf<-sim.vam(Time & Type ~ (ARAInf(0.5)|Weibull(1,3)))
simData<-simulate(simARAInf,5)

## plot1 ##
simulate(simARAInf,5)
plot(simARAInf,"v",col="blue")