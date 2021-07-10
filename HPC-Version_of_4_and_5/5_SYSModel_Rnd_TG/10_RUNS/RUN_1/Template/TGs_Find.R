## Last Modification 5/31/2021
## Find Target Elements of a Gene

TGs_Find<-function(GI){
  IndxTGs<-which(!is.na(match(Influ_net$GeneA,GI))) ## Ceck if we have this Gene in the network
  if(length(IndxTGs)){
    Gen_nam<-as.character(Influ_net$GeneB[IndxTGs])
    Conn_Types<-as.character(Influ_net$Connection[IndxTGs])
    TG_All<-data.frame("Elements"=Gen_nam,"Conn_Types"=Conn_Types)
    
    ####################### ALL Possible Connections
    ###################### To Family
    Family_Idx=which(grepl("(family)",as.character(TG_All$Elements), fixed = TRUE))
    
    if(length(Family_Idx)){  ### if there is at least one family element
      Curr_Family=TG_All[Family_Idx,]
      TG_All=TG_All[-Family_Idx,]
      if(length(which(Connection_Flag %in% Activa_Type))==0){  ### if there was not seen any -a 
        if(length(which(tail(Connection_Flag,1) %in% Trans_Types))==0){
          if(length(Curr_Family$Elements)){
            Family_Ele=Curr_Family
            colnames(Family_Ele)<-c("Complex_Elements","Con_Type")
            Curr_Family=NULL
          }else{
            Family_Ele=NULL
          }
        }else{ ##if Connection to this family has -t the next one can not be -t
          Trans_IDX=which(!is.na(match(Curr_Family$Conn_Types,Trans_Types))) ## remove family with -t connection
          if(length(Trans_IDX)){
            Curr_Family=Curr_Family[-Trans_IDX,]
          }
          if(length(Curr_Family$Conn_Types)){
            Family_Ele=Curr_Family
            colnames(Family_Ele)<-c("Complex_Elements","Con_Type")
            Curr_Family=NULL
          }else{
            Family_Ele=NULL
          } 
        }
      }else{ ### if there was any -a in the path to this element
        Activated_IDX=which(!is.na(match(Curr_Family$Conn_Types,Activa_Type)))
        if(length(Activated_IDX)){
          Curr_Family=Curr_Family[-Activated_IDX,]
        }
        if(length(which(tail(Connection_Flag,1) %in% Trans_Types))==0){
          if(length(Curr_Family$Elements)){
            Family_Ele=Curr_Family
            colnames(Family_Ele)<-c("Complex_Elements","Con_Type")
            Curr_Family=NULL
          }else{
            Family_Ele=NULL
          }
        }else{ ##if Connection to this family has -t the next one can not be -t
          Trans_IDX=which(!is.na(match(Curr_Family$Conn_Types,Trans_Types))) ## remove family with -t connection
          if(length(Trans_IDX)){
            Curr_Family=Curr_Family[-Trans_IDX,]
          }
          if(length(Curr_Family$Conn_Types)){
            Family_Ele=Curr_Family
            colnames(Family_Ele)<-c("Complex_Elements","Con_Type")
            Curr_Family=NULL
          }else{
            Family_Ele=NULL
          } 
        }
      }
    }else{
      Family_Ele=NULL
    }
 
    ###################### To Gene
    Genes_Idx=which(!is.na(match(as.character(TG_All$Elements),Exist_TG)))
    if((length(Genes_Idx)!=0)){
      Curr_Genes=TG_All[Genes_Idx,]
      TG_All=TG_All[-Genes_Idx,]
      if(length(which(Connection_Flag %in% Activa_Type))==0){
        Gen_Idx=which((!is.na(match(as.character(Curr_Genes$Conn_Types),Trans_Types))))
        if(length(Gen_Idx)){ ##if -t exist
          Gene_Ele=Curr_Genes[Gen_Idx,]
          colnames(Gene_Ele)<-c("Gene_Elements","Con_Type")
          Curr_Genes=Curr_Genes[-Gen_Idx,]
          if(length(Curr_Genes$Elements)){
            Connected_Genes=Curr_Genes
            colnames(Connected_Genes)<-c("Complex_Elements","Con_Type")
            Curr_Genes=NULL
          }else{
            Connected_Genes=NULL
          }
        }else{
          if(length(Curr_Genes$Elements)){
            Connected_Genes=Curr_Genes
            colnames(Connected_Genes)<-c("Complex_Elements","Con_Type")
          }else{
            Connected_Genes=NULL
          }
          Curr_Genes=NULL
          Gene_Ele=NULL
        }
      }else{
        Activated_IDX=which(!is.na(match(Curr_Genes$Conn_Types,Activa_Type)))
        if(length(Activated_IDX)){
          Curr_Genes=Curr_Genes[-Activated_IDX,]
          Connected_Genes=NULL
        }
        ### Find the target elements that are Gene and their connection is -t
        Gen_Idx=which((!is.na(match(Curr_Genes$Conn_Types,Trans_Types))))
        if(length(Gen_Idx)){ ##if -t exist
          Gene_Ele=Curr_Genes[Gen_Idx,]
          colnames(Gene_Ele)<-c("Gene_Elements","Con_Type")
          Curr_Genes=Curr_Genes[-Gen_Idx,]
          if(length(Curr_Genes$Elements)){
            Connected_Genes=Curr_Genes
            colnames(Connected_Genes)<-c("Complex_Elements","Con_Type")
            Curr_Genes=NULL
          }else{
            Connected_Genes=NULL
          }
        }else{
          if(length(Curr_Genes$Elements)){
            Connected_Genes=Curr_Genes
            colnames(Connected_Genes)<-c("Complex_Elements","Con_Type")
          }else{
            Connected_Genes=NULL
          }
          Curr_Genes=NULL
          Gene_Ele=NULL
        }
      }
    }else{
      Gene_Ele=NULL
      Connected_Genes=NULL
    }
    ####################### To a Complex, abstract and other elements which are not gene of complex or family or abstract. However they are in the network.
    Complex_Idx=which(grepl("(complex)",as.character(TG_All$Elements), fixed = TRUE)|grepl("(abstract)",as.character(TG_All$Elements), fixed = TRUE))
    if(length(Complex_Idx)){
      Curr_Complx=TG_All[Complex_Idx,]
      TG_All=TG_All[-Complex_Idx,]
      if(length(which(Connection_Flag%in% Activa_Type))==0){
        Complex_Elem=Curr_Complx
        colnames(Complex_Elem)<-c("Complex_Elements","Con_Type")
      }else{
        Activated_IDX=which(!is.na(match(Curr_Complx$Conn_Types,Activa_Type)))
        if(length(Activated_IDX)){
          Curr_Complx=Curr_Complx[-Activated_IDX,]
        }
        Complex_Rem_Idx=length(Curr_Complx$Elements)
        if(Complex_Rem_Idx){
          Complex_Elem=Curr_Complx
          colnames(Complex_Elem)<-c("Complex_Elements","Con_Type")
          Curr_Complx=Curr_Complx[-Complex_Rem_Idx,]
        }else{
          Complex_Elem=NULL
        }
      }
    }else{
      Complex_Elem=NULL
    }
    ############### Other Elements
    if(length(TG_All$Elements)){
      Curr_Other_Ele=TG_All
      if(length(which(Connection_Flag %in% Activa_Type))==0){
        Other_ELe_Activate=which(Curr_Other_Ele$Conn_Types %in% Activa_Type)
        if(length(Other_ELe_Activate)!=0){
          Other_Elemnts=Curr_Other_Ele[Other_ELe_Activate,]
          colnames(Other_Elemnts)<-c("Complex_Elements","Con_Type")
          TG_All=NULL
        }else{
          TG_All=NULL
          Other_Elemnts=NULL
        }
      }else{
        TG_All=NULL
        Other_Elemnts=NULL
      }
    }else{
      Other_Elemnts=NULL
    }
    
    ##################
    Other_ELEMS_List<-list(Complex_Elem,Family_Ele,Connected_Genes,Other_Elemnts)
    Complex_Ele = do.call(rbind, Other_ELEMS_List)
    ELEMS<-list(Gene_Ele,Complex_Ele)
    return(ELEMS)
  }else{
    #cat("!!!!Current Gene of Interest does not connected to the other Gene in our Influence Network, GeneNAme:",GI,"\n")
    return(NULL)
  }
}