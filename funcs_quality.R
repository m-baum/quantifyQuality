# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Auxiliary functions for                                                     #
# Quantifying the Quality of Configurational Causal Models                    #
#                                                                             #
# [anonymized]                                                                #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


candidates <- function(x){
  ## Preparation
  x <- noblanks(x)
  asfsX <- unlist(extract_asf(x))
  asfsX <- lapply(asfsX, function(z) paste0(stdCond(lhs(z)),"<->", rhs(z)))

  if(any(unlist(lapply(str_split(lhs(asfsX),"[\\+\\*]"), function(l) lapply(l, nchar)))!=1))
    stop("Please use only single letters as variable/factor names.")
  
  if(length(asfsX)==1){
    liCan <- unique(unlist(str_split(lhs(unlist(asfsX)),"[*,+]")))
    lit.cand <- split(liCan,rhs(asfsX))
    con.cand <- str_split(lhs(asfsX), "\\+")
    names(con.cand) <- rhs(asfsX)
    dis.cand <- split(lhs(asfsX),rhs(asfsX))
    vert.cand <- split(lapply(liCan, function(l) append(l,rhs(asfsX))),rhs(asfsX))
    # names(vert.cand[[1]]) <- rhs(asfsX)
   for(i in 1:length(vert.cand)){
     names(vert.cand[[i]]) <- rep(names(vert.cand[i]), length(vert.cand[[i]]))
    
    }
    return(list(literal.candidates=lit.cand, conjunction.candidates=con.cand, disjunction.candidates=dis.cand, sequential.candidates=vert.cand))
  }
  
  ## Literal Candidates (for single letter factors)
  directC <- unlist(lapply(split(asfsX, rhs(asfsX)), function(t) str_extract_all(lhs(t),"[A-Za-z]")), recursive = F)
  eff <- names(directC)
  
  directC_NEG <- lapply(directC, switch.case)
  names(directC_NEG) <- switch.case(eff)
  relations <- c(directC,directC_NEG)
  
  ## Chain Search (for single letter factors)
  chains <- as.list(eff)
  end<-FALSE
  while(end==FALSE){
    for(i in 1:length(chains)){
      lastElement <- str_sub(chains[[i]], nchar(chains[[i]]), nchar(chains[[i]]))
      if(lastElement%in%names(relations)){
        chain <- lapply(unname(unlist(relations[which(names(relations)==lastElement)])), function(l) paste0(chains[[i]],l))
        chains[[i]] <- unlist(chain)
      }else if(lastElement != "-"){
          chains[[i]] <- paste0(chains[[i]],"-")}
      }
    chains <- as.list(unlist(chains))
    if(all(lapply(chains, function(s) str_sub(s, nchar(s), nchar(s)))=="-")==TRUE)
      end <- TRUE
    }
  chains <- unlist(lapply(chains, function(s) str_remove(s, "-"))) #to be read backwards, i.e. element on the very left of a string is the last element of the chain
  
  relLit.cand <- split(chains, str_sub(chains,1,1))
  relLit.cand <- lapply(relLit.cand, function(t) str_remove(t, str_sub(t,1,1)))
  relLit.cand <- lapply(relLit.cand, function(r) unique(unlist(strsplit(r,split=""))))
  
  # Disjunction Candidates (disjunctive connection of all conjunction candidates):
  #cheap version: disj.cand <- lapply(conjunction.cand, function(c) paste0(c, collapse = "+"))
  
  ## Disjunction Candidates
  asfsX_con <- lapply(split(unlist(asfsX),rhs(unlist(asfsX))), lhs)
  asfsX_con_NEG <- lapply(asfsX_con, negate3)
  names(asfsX_con_NEG) <- switch.case(names(asfsX_con))
  
  conjunc <- lapply(asfsX_con, function(k) unlist(str_split(k, "\\+"), recursive=F))
  conjunc_NEG <- lapply(asfsX_con_NEG, function(k) unlist(str_split(k, "\\+"), recursive=F))
  conjunc_all <- c(conjunc,conjunc_NEG)
  
  pot_sub <- names(conjunc)
  Disj.cand <- lapply(conjunc, list)
  
  for(k in 1:length(Disj.cand)){
    if(any(str_detect(paste0(unlist(Disj.cand[[k]]),collapse="*"), names(conjunc_all)))){
      dis.can.temp <- Disj.cand[[k]]
      end<-FALSE
      while(end==FALSE){
        for(i in 1:length(dis.can.temp)){
          if(any(str_detect(paste0(unlist(dis.can.temp[[i]]),collapse="*"),names(conjunc_all)))){
            subs <- lapply(pot_sub, function(s) 
              if(any(str_detect(dis.can.temp[[i]],s))||
                 any(str_detect(dis.can.temp[[i]],switch.case(s)))){
                unlist(lapply(dis.can.temp[[i]], function(k) if(str_detect(k,s)){
                  str_replace(k, s, unlist(conjunc_all[which(names(conjunc_all)==s)]))
                }else if(str_detect(k,switch.case(s))){
                  str_replace(k, switch.case(s), unlist(conjunc_all[which(names(conjunc_all)==switch.case(s))]))
                }else{k}))
              })
            subs <- subs[-which(sapply(subs, is.null))]
            dis.can.temp[[i]] <- lapply(subs, function(y) unlist(lapply(y, function(z) paste0(unique(unlist(str_split(z, "\\*"))),collapse = "*"))))
          }else{
            dis.can.temp[[i]] <- as.list(is.null)
          }
        }
        dis.can.temp <- unlist(dis.can.temp, recursive = F)
        dis.can.temp <- dis.can.temp[!sapply(dis.can.temp,is.null)]
        #dis.can.temp <- lapply(dis.can.temp, unique)
        Disj.cand[[k]] <- append(Disj.cand[[k]],dis.can.temp)
        
        if(length(dis.can.temp)==0) end <- TRUE
      }
    }
  }
  
  dnfs <- lapply(Disj.cand, function(t) lapply(t, function(l) stdCond(paste0(l, collapse = "+"))))
  disj.cand <- lapply(dnfs, function(s) unique(stdCond(unlist(s))))
  
  ## Conjunction Candidates
  conjunction.cand <- lapply(Disj.cand, function(m) unique(stdCond(unlist(m))))
  
  # Sequential Candidates
  # (note: Indirect relations in paths are to be checked in relevancies() with corresponding asfs. 
  # If indirect relevancies do not exist, i.e. if there is no such asf, the entire path is deleted.)
  vertical.cand <- unlist(lapply(chains, stri_reverse))
  vertical.cand <- split(vertical.cand, substr(vertical.cand, nchar(vertical.cand), nchar(vertical.cand)))
  vertical.cand <- lapply(vertical.cand, function(l) lapply(l, function(j) unlist(strsplit(j,split=""))))
  
  lapply(vertical.cand, function(x) unlist(x, recursive = F))
  
  for(i in 1:length(vertical.cand)){
     names(vertical.cand[[i]]) <- rep(names(vertical.cand[i]), length(vertical.cand[[i]]))
    }
  
  # Result
  list(literal.candidates=relLit.cand, conjunction.candidates=conjunction.cand, disjunction.candidates=disj.cand, sequential.candidates=vertical.cand)
}



relevancies <- function(x){
  # ### Preparation
  if(any(class(x)=="character")){
    csfs <- x
    candidates.x <- lapply(x, candidates)
    x <- noblanks(x)
    asfs <- extract_asf(x)
    outcomes <- lapply(asfs, rhs)
    data <- lapply(x, function(y) selectCases(y))
    
    # using CNA 
    max.comp <- max(mapply(function(x,y) ncol(x)-length(y),x=data,y=outcomes))
    ana <- mapply(function(x,y) suppressMessages(cna(x,outcome=y, rm.dup.factors = F,
                                                     maxstep = c(max.comp,2*max.comp,3*(max.comp)+5))), x= data,y=outcomes, SIMPLIFY = F)  
    
    asfs.ana <- lapply(ana, function(x) asf(x)$condition)

  }
  if(any(class(x)=="cna")){
    csfs <- suppressWarnings(csf(x)$condition)
    if(length(csfs)==0){csfs <- suppressWarnings(asf(x)$condition)}
    candidates.x <- lapply(csfs, candidates)
    outcomes <- lapply(csfs, function(x) extract_asf(x))
    outcomes <- lapply(outcomes, function(x) lapply(x, function(y) rhs(y)))
    asfs <- asf(x)
    asfs.ana <- vector("list", length(csfs))
    for(i in 1:length(asfs.ana)){
      asfs.ana[[i]] <-  asfs[asfs$outcome%in%unlist(outcomes[[i]]),]$condition
      
    }
    # asfs.ana <- lapply(asfs.ana, function(s) s <- asfs)
  }
  
  # Check if recovered asfs are submodels of candidate disjunctions
  relevant.asfs <- vector("list", length(asfs.ana))
  for(i in 1:length(asfs.ana)){
    
     recovered.asfs <- split(asfs.ana[[i]],rhs(asfs.ana[[i]]))
     cand.asfs <- candidates.x[[i]]$disjunction.candidates
     for(j in 1:length(cand.asfs)){
       cand.asfs[[j]] <-  paste0(cand.asfs[[j]],"<->",names(cand.asfs[j]))
     }
     
     score.j <- vector("list",length(recovered.asfs))
     for(j in 1:length(recovered.asfs)){
     score.k <- vector("list",length(recovered.asfs[[j]]))
        for(k in 1:length(recovered.asfs[[j]])){
       
        rr <- expand.grid(recovered.asfs[[j]][k],cand.asfs[[j]])
        score.k[[k]] <- any(apply(rr,1,function(x) is.submodel(x[1],x[2])))
        }
       score.j[[j]] <- unlist(score.k)
        }
       
    relevant.asfs[[i]] <-   mapply(function(x,y) x[y],x=recovered.asfs,y=score.j)
     
  }
  
  
  #### Extract relevant literals, conjunctions, disjunctions, and sequences from relevant.asfs
    relevant.Literals.x <- vector("list", length(relevant.asfs))
  relevant.Conjuncts.x <- vector("list", length(relevant.asfs))
  relevant.Disjunctions.x <- vector("list", length(relevant.asfs))
  relevant.Verticals.x <- vector("list", length(relevant.asfs))
  for(i in 1:length(relevant.asfs)){
    conj.ana <- lapply(relevant.asfs[[i]], function(x) unlist(str_split(lhs(x),"[+]")))
    lit.ana <- lapply(conj.ana, function(x) unlist(str_split(x,"[*]")))
    relevant.Literals.x[[i]] <- lapply(lit.ana, function(x) unique(x))
    
    relevant.Conjuncts.x[[i]] <- lapply(conj.ana, function(x) unique(x))
    
    # disjunctions
    relevant.Disjunctions.x[[i]] <- lapply(relevant.asfs[[i]], function(x) lhs(x))

    # identify relevant sequences 
    select3 <-      candidates.x[[i]]$sequential.candidates
    r.Vert <- vector("list", length(select3))
    names(r.Vert) <- names(candidates.x[[i]]$sequential.candidates)
    for(j in 1:length(select3)){
      aux.vert <- vector("list", length(select3[[j]]))
      for(k in 1:length(select3[[j]])){
        test <-   select3[[j]][[k]]
        aux.vert[[k]]  <-  all(test[-length(test)] %in% unique(unlist(relevant.Literals.x[[i]][which(names(relevant.Literals.x[[i]])%in% test[length(test)])])))
      }
      r.Vert[[j]] <- select3[[j]][unlist(aux.vert)]
      r.Vert[[j]] <- r.Vert[[j]][!duplicated(r.Vert[[j]])]
      
    }
    
    relevant.Verticals.x[[i]] <- r.Vert
    
    relevant.Literals.x[[i]] <-  relevant.Literals.x[[i]][unlist(outcomes[[i]])]
    relevant.Conjuncts.x[[i]] <-  relevant.Conjuncts.x[[i]][unlist(outcomes[[i]])]
    relevant.Disjunctions.x[[i]] <-  relevant.Disjunctions.x[[i]][unlist(outcomes[[i]])]
    relevant.Verticals.x[[i]] <-  relevant.Verticals.x[[i]][unlist(outcomes[[i]])]
  }
 
  
  ########
  # Prepare output
  if(any(class(x)=="character")){
    for(i in 1:length(x)){
      attr(relevant.Literals.x[[i]], "model1") <- x[i]
      attr(relevant.Conjuncts.x[[i]], "model2") <- x[i]
      attr(relevant.Disjunctions.x[[i]], "model3") <- x[i]
      attr(relevant.Verticals.x[[i]], "model4") <- x[i]
    }
    out <- list(relevant.Literals=relevant.Literals.x , relevant.Conjuncts=relevant.Conjuncts.x, 
                relevant.Disjunctions=relevant.Disjunctions.x, relevant.Sequences=relevant.Verticals.x, ana=ana,candidates=candidates.x)
    out <- structure(out, class = "rel")
  }
  if(any(class(x)=="cna")){  
    for(i in 1:length(csfs)){
      attr(relevant.Literals.x[[i]], "model1") <- csfs[i]
      attr(relevant.Conjuncts.x[[i]], "model2") <- csfs[i]
      attr(relevant.Disjunctions.x[[i]], "model3") <- csfs[i]
      attr(relevant.Verticals.x[[i]], "model4") <- csfs[i]
    }
    out <- list(relevant.Literals=relevant.Literals.x , relevant.Conjuncts=relevant.Conjuncts.x, 
                relevant.Disjunctions=relevant.Disjunctions.x, relevant.Sequences=relevant.Verticals.x, csfs=csfs, candidates=candidates.x)
    out <- structure(out, class = "rel")
  }
  out
}


print.rel <- function(x, ...){
  print(x[1:4])
}




quality <- function(x, GT, beta=1, agregate="mean1"){
  
  #===== PREPARATION 
  #### Identify relevant literals in x, if x is not a cna object
  if(any(class(x)=="character")){
  y <- relevancies(x)
  # relevant.Literals.M <- y$relevant.Literals
  ana.M <- y$ana
  # candidates.M <- y$candidates
  }
  ###### # Identify relevant literals in x, if x is a cna object
  if(any(class(x)=="cna")){
  y <- relevancies(x)
  # relevant.Literals.M <- y$relevant.Literals
  csfs.M <- y$csfs
  # candidates.M <- y$candidates
  }
  
  ########
  
  ##### Identify relevant literals in GT
  yy <- relevancies(GT) 
  relevant.Literals.G <- yy$relevant.Literals[[1]] 

  
  ###### LITERAL QUALITY

  # ====== Literal correctness
  relevant.Literals.M <- y$relevant.Literals
  litCor <- vector("list", length(relevant.Literals.M))
  for(i in 1:length(relevant.Literals.M)){
    # cat(i, "\n")
    aux <- vector("list", length(relevant.Literals.M[[i]]))
    for(j in 1:length(relevant.Literals.M[[i]])){
      fac <- names(relevant.Literals.M[[i]])[j]
      compare <- relevant.Literals.G[which(names(relevant.Literals.G) %in% fac)]
      if(length(compare)>0){
        inter <- intersect(relevant.Literals.M[[i]][[j]], unlist(compare))
        score <- length(inter)/length(relevant.Literals.M[[i]][[j]])
        weight <- length(relevant.Literals.M[[i]][[j]])/length(unlist(relevant.Literals.M[[i]]))
        # rel.lit <- I(list(relevant.Literals.M[[i]][[j]]))
        aux[[j]] <- data.frame(Cor.score=score,weight=weight)# ,rel.literals=rel.lit)
      }else{
        score <-  0
        weight <- length(relevant.Literals.M[[i]][[j]])/length(unlist(relevant.Literals.M[[i]]))
        # rel.lit <- I(list(relevant.Literals.M[[i]][[j]]))
        aux[[j]] <-data.frame(Cor.score=score,weight=weight)
      }
    }
    litCor[[i]] <- aux
    names(litCor[[i]]) <- names(relevant.Literals.M[[i]])
    litCor[[i]] <- do.call(rbind, litCor[[i]])
    # if(agregate=="harmonic2"){
    #    litCor[[i]]$weighted.mean <- sum(litCor[[i]]$weight) / sum(litCor[[i]]$weight / litCor[[i]]$Cor.score)
    # }else{
    litCor[[i]]$weighted.mean <- weighted.mean(litCor[[i]]$Cor.score,litCor[[i]]$weight)
    # }
    litCor[[i]] <- list(litCor[[i]])
    names(litCor[[i]]) <- "Literal.Correctness"
     attr(litCor[[i]], "relevant.Literals") <- relevant.Literals.M[[i]]
  }
  
  
  # ====== Literal completeness
  litCom <- vector("list", length(relevant.Literals.M))
  for(i in 1:length(relevant.Literals.M)){
    # cat(i, "\n")
    aux <- vector("list", length(relevant.Literals.G))
    for(j in 1:length(relevant.Literals.G)){
      fac <- names(relevant.Literals.G[j])
      compare <- relevant.Literals.M[[i]][which(names(relevant.Literals.M[[i]]) %in% fac)]
      if(length(compare)>0){
        inter <- intersect(relevant.Literals.G[[j]], unlist(compare))
        score <- length(inter)/length(relevant.Literals.G[[j]])
        weight <- length(relevant.Literals.G[[j]])/length(unlist(relevant.Literals.G))
        # rel.lit <- I(list(relevant.Literals.M[[i]][[j]]))
        aux[[j]] <- data.frame(Com.score=score,weight=weight)# ,rel.literals=rel.lit)
      }else{
        score <-  0
        weight <- length(relevant.Literals.G[[j]])/length(unlist(relevant.Literals.G))
        # rel.lit <- I(list(relevant.Literals.M[[i]][[j]]))
        aux[[j]] <-data.frame(Com.score=score,weight=weight)
      }
    }
    litCom[[i]] <- aux
    names(litCom[[i]]) <- names(relevant.Literals.G)
    litCom[[i]] <- do.call(rbind, litCom[[i]])
    litCom[[i]]$weighted.mean <- weighted.mean(litCom[[i]]$Com.score,litCom[[i]]$weight)
    litCom[[i]] <- list(litCom[[i]])
     names(litCom[[i]]) <- "Literal.Completeness"
     attr(litCom[[i]], "relevant.Literals") <- relevant.Literals.G
  }
  
   # ==== FBeta Aggregation
  
  FBeta.Lit <- vector("list", length(relevant.Literals.M))
  for(i in 1:length(FBeta.Lit)){
    FBeta <- (1+beta^2)*(litCor[[i]][[1]]$weighted.mean[1] * litCom[[i]][[1]]$weighted.mean[1] )/((beta^2*litCor[[i]][[1]]$weighted.mean[1])+litCom[[i]][[1]]$weighted.mean[1]) 
    if(is.na(FBeta)){FBeta <- 0}
    FBeta.Lit[[i]] <-  list(litCor[[i]],litCom[[i]],FBeta=FBeta)
    attr(FBeta.Lit[[i]], "model") <- attr(relevant.Literals.M[[i]], "model1")
    attr(FBeta.Lit[[i]], "groundTruth") <- GT
}
  
  
  
  ####### CONJUNCTIVE  QUALITY
  
  # Relevant conjuncts of M
   relevant.Conjuncts.M <- y$relevant.Conjuncts
  # Relevant conjuncts of GT
   relevant.Conjuncts.G <- yy$relevant.Conjuncts[[1]]
   
   
   # ======= Conjunctive Correctness
  ConjCor <- vector("list", length(relevant.Conjuncts.M))
  for(i in 1:length(relevant.Conjuncts.M)){
    # cat(i, "\n")
    aux <- vector("list", length(relevant.Conjuncts.M[[i]]))
    for(j in 1:length(relevant.Conjuncts.M[[i]])){
      fac <- names(relevant.Conjuncts.M[[i]])[j]
      compare <- relevant.Conjuncts.G[which(names(relevant.Conjuncts.G) %in% fac)]
      if(length(compare)>0){
        score <- vector("list", length(relevant.Conjuncts.M[[i]][[j]]))
        for(k in 1:length(relevant.Conjuncts.M[[i]][[j]])){
        expan <- expand.grid(relevant.Conjuncts.M[[i]][[j]][k], unlist(compare),stringsAsFactors = F)
        aux3 <- vector("list", nrow(expan))
        for(g in 1:nrow(expan)){
          # aux3[[g]] <- is.submodel(paste0(expan[g,1],"<->X"),paste0(expan[g,2],"<->X"))
          aux3[[g]] <- asf.overlap4(c(paste0(expan[g,1],"<->X"),paste0(expan[g,2],"<->X")))$max.intersect
        }
        max.con <- unlist(aux3)
        q <- str_extract_all(max.con, "[*]") 
        q <- unlist(lapply(q, length))
        # if(length(unlist(str_extract_all(relevant.Conjuncts.M[[i]][[j]][k],"[*]")))==0){aux2 <- NA}else{
        if(suppressWarnings(max(q)>0)){aux2 <- max.con[which(q==max(q)&nchar(max.con)==max(nchar(max.con)))][1]}else{
            if(suppressWarnings(length(max.con)>0)){aux2 <- max.con[which(nchar(max.con)==max(nchar(max.con)))][1]}else{
             aux2  <- 0 
            }
        }
        # }

          if(!is.numeric(aux2)){
        score[[k]]  <- getComplexity(aux2)/getComplexity(relevant.Conjuncts.M[[i]][[j]][k])
        }else{score[[k]]  <- 0}
        # }
        }
        
        w_weights <- getComplexity(relevant.Conjuncts.M[[i]][[j]])/sum(getComplexity(relevant.Conjuncts.M[[i]][[j]]))
        
        score.overall <- weighted.mean(unlist(score),w_weights,na.rm=T)

        weight <- sum(getComplexity(relevant.Conjuncts.M[[i]][[j]]))/sum(getComplexity(unlist(relevant.Conjuncts.M[[i]])))
        # rel.con <- I(list(relevant.Conjuncts.M[[i]][[j]]))
        aux[[j]] <- data.frame(Cor.score=score.overall,weight=weight)# ,relevant.conjunctions=rel.con)
      }else{
        score.overall <-  0
         weight <- sum(getComplexity(relevant.Conjuncts.M[[i]][[j]]))/sum(getComplexity(unlist(relevant.Conjuncts.M[[i]])))
        aux[[j]] <-data.frame(Cor.score=score.overall,weight=weight)# ,relevant.conjunctions=rel.con)
      }
    }
    ConjCor[[i]] <- aux
    names(ConjCor[[i]]) <- names(relevant.Conjuncts.M[[i]])
    ConjCor[[i]] <- do.call(rbind, ConjCor[[i]])
    ConjCor[[i]]$weighted.mean <- weighted.mean(ConjCor[[i]]$Cor.score,ConjCor[[i]]$weight,na.rm = T)
    ConjCor[[i]] <- list(ConjCor[[i]])
    names(ConjCor[[i]]) <- "Conjunction.Correctness"
     attr(ConjCor[[i]], "relevant.Conjuncts") <- relevant.Conjuncts.M[[i]]
  }
    
  
   # ======= Conjunctive Completeness
  ConjCom <- vector("list", length(relevant.Conjuncts.M))
    for(i in 1:length(relevant.Conjuncts.M)){
    # cat(i, "\n")
    aux <- vector("list", length(relevant.Conjuncts.G))
    for(j in 1:length(relevant.Conjuncts.G)){
      fac <- names(relevant.Conjuncts.G)[j]
      compare <- unlist(relevant.Conjuncts.M[[i]][which(names(relevant.Conjuncts.M[[i]]) %in% fac)])
      if(length(compare)>0){
        score <- vector("list", length(compare))
        for(k in 1:length(relevant.Conjuncts.G[[j]])){
        expan <- expand.grid(relevant.Conjuncts.G[[j]][k], unlist(compare) ,stringsAsFactors = F)
        aux3 <- vector("list", nrow(expan))
        for(g in 1:nrow(expan)){
          # aux3[[g]] <- is.submodel(paste0(expan[g,1],"<->X"),paste0(expan[g,2],"<->X"))
          aux3[[g]] <- asf.overlap4(c(paste0(expan[g,1],"<->X"),paste0(expan[g,2],"<->X")))$max.intersect
        }
        max.con <- unlist(aux3)
        q <- str_extract_all(max.con, "[*]") 
        q <- unlist(lapply(q, length))
        if(suppressWarnings(max(q)>0)){aux2 <- max.con[which(q==max(q)&nchar(max.con)==max(nchar(max.con)))][1]}else{
            if(suppressWarnings(length(max.con)>0)){aux2 <- max.con[which(nchar(max.con)==max(nchar(max.con)))][1]}else{
             aux2  <- 0 
            }
        }
        

          if(!is.numeric(aux2)){
        score[[k]]  <- getComplexity(aux2)/getComplexity(relevant.Conjuncts.G[[j]][k])
        }else{score[[k]]  <- 0}
        # }
        }
        w_weights <- getComplexity(relevant.Conjuncts.G[[j]])/sum(getComplexity(relevant.Conjuncts.G[[j]]))
        
        score.overall <- weighted.mean(unlist(score),w_weights,na.rm=T)
        
        # score.overall <- mean(unlist(score),na.rm = T)

        weight <- sum(getComplexity(unlist(relevant.Conjuncts.G[j])))/sum(getComplexity(unlist(relevant.Conjuncts.G)))
        # rel.con <- I(list(relevant.Conjuncts.M[[i]][[j]]))
        aux[[j]] <- data.frame(Com.score=score.overall,weight=weight)# ,relevant.conjunctions=rel.con)
      }else{
        score.overall <-  0
        weight <- sum(getComplexity(unlist(relevant.Conjuncts.G[j])))/sum(getComplexity(unlist(relevant.Conjuncts.G)))
        aux[[j]] <-data.frame(Com.score=score.overall,weight=weight)# ,relevant.conjunctions=rel.con)
      }
    }
    ConjCom[[i]] <- aux
    names(ConjCom[[i]]) <- names(relevant.Conjuncts.G)
    ConjCom[[i]] <- do.call(rbind, ConjCom[[i]])
    ConjCom[[i]]$weighted.mean <- weighted.mean(ConjCom[[i]]$Com.score,ConjCom[[i]]$weight,na.rm = T)
    ConjCom[[i]] <- list(ConjCom[[i]])
    names(ConjCom[[i]]) <- "Conjunction.Completeness"
     attr(ConjCom[[i]], "relevant.Conjuncts") <- relevant.Conjuncts.G
  }
  
  
 
  
  
  # ====== FBeta Aggregation
  
    FBeta.Conj <- vector("list", length(relevant.Conjuncts.M))
    for(i in 1:length(FBeta.Conj)){
      FBeta <- (1+beta^2)*(ConjCor[[i]][[1]]$weighted.mean[1] * ConjCom[[i]][[1]]$weighted.mean[1] )/((beta^2*ConjCor[[i]][[1]]$weighted.mean[1])+ConjCom[[i]][[1]]$weighted.mean[1]) 
      if(is.na(FBeta)){FBeta <- 0}
      FBeta.Conj[[i]] <-  list(ConjCor[[i]],ConjCom[[i]],FBeta=FBeta)
      attr(FBeta.Conj[[i]], "model") <- attr(relevant.Conjuncts.M[[i]], "model2")
      attr(FBeta.Conj[[i]], "groundTruth") <- GT
      }
  
  
  
  
  ####### DISJUNCTIVE QUALITY
  
  
  # Relevant disjuncts of M
   relevant.Disjunctions.M <- y$relevant.Disjunctions
  # Relevant conjuncts of GT
   relevant.Disjunctions.G <- yy$relevant.Disjunctions[[1]]

  # ======= Disjunctive correctness
   DisjCor <- vector("list", length(relevant.Disjunctions.M))
  for(i in 1:length(relevant.Disjunctions.M)){
    # cat(i, "\n")
    aux <- vector("list", length(relevant.Disjunctions.M[[i]]))
    for(j in 1:length(relevant.Disjunctions.M[[i]])){
      fac <- names(relevant.Disjunctions.M[[i]])[j]
      compare <- relevant.Disjunctions.G[which(names(relevant.Disjunctions.G) %in% fac)]
      if(length(unlist(compare))>0){
        # aux2 <- vector("list", length(relevant.Disjunctions.M[[i]][[j]]))
        score <- vector("list", length(relevant.Disjunctions.M[[i]][[j]]))

        for(k in 1:length(relevant.Disjunctions.M[[i]][[j]])){
        expan <- expand.grid(relevant.Disjunctions.M[[i]][[j]][k], unlist(compare),stringsAsFactors = F)
        aux3 <- vector("list", nrow(expan))
        for(g in 1:nrow(expan)){
          aux3[[g]] <- asf.overlap4(c(paste0(expan[g,1],"<->X"),paste0(expan[g,2],"<->X")))$max.intersect
        }
        max.dis <- unlist(aux3)
        q <- str_extract_all(max.dis, "[+]") 
        q <- unlist(lapply(q, length))
        # if(length(unlist(str_extract_all(relevant.Disjunctions.M[[i]][[j]][k],"[+]")))==0){aux2 <- NA}else{
          if(suppressWarnings(max(q)>0)){aux2 <- max.dis[which(q==max(q)&nchar(max.dis)==max(nchar(max.dis)))][1]}else{
          if(suppressWarnings(max(q)==0)){
        # aux4 <- vector("list", nrow(expan))
        # for(g in 1:nrow(expan)){
        #   aux4[[g]] <- asf.overlap4(c(paste0(expan[g,1],"<->X"),paste0(expan[g,2],"<->X")))$max.intersect
        # }
        aux2  <- max.dis[which(q==max(q)&nchar(max.dis)==max(nchar(max.dis)))][1] 
          }else{
            aux2  <- 0}
        }
        # }
        

          if(!is.numeric(aux2)){
        score[[k]]  <- getComplexity(aux2)/sum(getComplexity(relevant.Disjunctions.M[[i]][[j]][k]))
        }else{score[[k]]  <- 0}
        # }
        }
        
        w_weights <- getComplexity(relevant.Disjunctions.M[[i]][[j]])/sum(getComplexity(relevant.Disjunctions.M[[i]][[j]]))
        
        score.overall <- weighted.mean(unlist(score),w_weights,na.rm=T)
        
        # score.overall <- mean(unlist(score),na.rm = T)
        weight <- sum(getComplexity(relevant.Disjunctions.M[[i]][[j]]))/sum(getComplexity(unlist(relevant.Disjunctions.M[[i]])))
        # rel.con <- I(list(relevant.Conjuncts.M[[i]][[j]]))
        aux[[j]] <- data.frame(Cor.score=score.overall,weight=weight)# ,relevant.conjunctions=rel.con)
      }else{
        score.overall <-  0
        weight <- sum(getComplexity(relevant.Disjunctions.M[[i]][[j]]))/sum(getComplexity(unlist(relevant.Disjunctions.M[[i]])))
        aux[[j]] <-data.frame(Cor.score=score.overall,weight=weight)# ,relevant.conjunctions=rel.con)
      }
    }
    DisjCor[[i]] <- aux
    names(DisjCor[[i]]) <- names(relevant.Disjunctions.M[[i]])
    DisjCor[[i]] <- do.call(rbind, DisjCor[[i]])
    DisjCor[[i]]$weighted.mean <- weighted.mean(DisjCor[[i]]$Cor.score,DisjCor[[i]]$weight,na.rm = T)
    DisjCor[[i]] <- list(DisjCor[[i]])
    names(DisjCor[[i]]) <- "Disjunction.Correctness"
    # attr(ConjCor[[i]], "model") <- attr(relevant.Conjuncts.M[[i]], "model")
    # attr(DisjCor[[i]], "groundTruth") <- GT
    attr(DisjCor[[i]], "relevant.Disjunctions") <- relevant.Disjunctions.M[[i]]
  }
  
  
  # ======= Disjunctive completeness
   DisjCom <- vector("list", length(relevant.Disjunctions.M))
  for(i in 1:length(relevant.Disjunctions.M)){
    # cat(i, "\n")
    aux <- vector("list", length(relevant.Disjunctions.G))
    for(j in 1:length(relevant.Disjunctions.G)){
      fac <- names(relevant.Disjunctions.G)[j]
      compare <- unlist(relevant.Disjunctions.M[[i]][which(names(relevant.Disjunctions.M[[i]]) %in% fac)])
      if(length(compare)>0){
        # aux2 <- vector("list", length(relevant.Disjunctions.M[[i]][[j]]))
        score <- vector("list", length(compare))

        for(k in 1:length(relevant.Disjunctions.G[[j]])){
        expan <- expand.grid(relevant.Disjunctions.G[[j]][k], unlist(compare),stringsAsFactors = F)
        aux3 <- vector("list", nrow(expan))
        for(g in 1:nrow(expan)){
          aux3[[g]] <- asf.overlap4(c(paste0(expan[g,1],"<->X"),paste0(expan[g,2],"<->X")))$max.intersect
        }
        max.dis <- unlist(aux3)
        q <- str_extract_all(max.dis, "[+]") 
        q <- unlist(lapply(q, length))
        if(suppressWarnings(max(q)>0)){aux2 <- max.dis[which(q==max(q)&nchar(max.dis)==max(nchar(max.dis)))][1]}else{
          if(suppressWarnings(max(q)==0)){
        aux2  <- max.dis[which(q==max(q)&nchar(max.dis)==max(nchar(max.dis)))][1]
          }else{aux2  <- 0}
        }

          if(!is.numeric(aux2)){
        score[[k]]  <- sum(getComplexity(aux2))/getComplexity(relevant.Disjunctions.G[[j]][k])
        }else{score[[k]]  <- 0}
        # }
        }
        
        w_weights <- getComplexity(relevant.Disjunctions.G[[j]])/sum(getComplexity(relevant.Disjunctions.G[[j]]))
        
        score.overall <- weighted.mean(unlist(score),w_weights,na.rm=T)
        
        # score.overall <- mean(unlist(score),na.rm = T)
        weight <- sum(getComplexity(unlist(relevant.Disjunctions.G[j])))/sum(getComplexity(unlist(relevant.Disjunctions.G)))
        aux[[j]] <- data.frame(Com.score=score.overall,weight=weight)# ,relevant.conjunctions=rel.con)
      }else{
        score.overall <-  0
        weight <- sum(getComplexity(unlist(relevant.Disjunctions.G[j])))/sum(getComplexity(unlist(relevant.Disjunctions.G)))
        # weight <- sum(getComplexity(unlist(compare)))/sum(getComplexity(unlist(relevant.Disjunctions.G)))
        aux[[j]] <-data.frame(Com.score=score.overall,weight=weight)# ,relevant.conjunctions=rel.con)
      }
    }
    DisjCom[[i]] <- aux
    names(DisjCom[[i]]) <- names(relevant.Disjunctions.G)
    DisjCom[[i]] <- do.call(rbind, DisjCom[[i]])
    DisjCom[[i]]$weighted.mean <- weighted.mean(DisjCom[[i]]$Com.score,DisjCom[[i]]$weight,na.rm = T)
    DisjCom[[i]] <- list(DisjCom[[i]])
    names(DisjCom[[i]]) <- "Disjunction.Completeness"
    # attr(ConjCor[[i]], "model") <- attr(relevant.Conjuncts.M[[i]], "model")
    # attr(DisjCor[[i]], "groundTruth") <- GT
    attr(DisjCom[[i]], "relevant.Disjunctions") <- relevant.Disjunctions.G
  }
  
  
  # ====== FBeta Aggregation
  
    FBeta.Disj <- vector("list", length(relevant.Disjunctions.M))
    for(i in 1:length(FBeta.Disj)){
      FBeta <- (1+beta^2)*(DisjCor[[i]][[1]]$weighted.mean[1] * DisjCom[[i]][[1]]$weighted.mean[1] )/
                  ((beta^2*DisjCor[[i]][[1]]$weighted.mean[1])+DisjCom[[i]][[1]]$weighted.mean[1]) 
      if(is.na(FBeta)){FBeta <- 0}
      FBeta.Disj[[i]] <-  list(DisjCor[[i]],DisjCom[[i]],FBeta=FBeta)
      attr(FBeta.Disj[[i]], "model") <- attr(relevant.Disjunctions.M[[i]], "model3")
      attr(FBeta.Disj[[i]], "groundTruth") <- GT
      }  
  
  
  
  # VERTICAL Quality
  relevant.Verticals.M <- y$relevant.Sequences
   # Relevant verticals of GT
  relevant.Verticals.G <- yy$relevant.Sequences[[1]]
  
  # ======= Vertical Correctness


   VertCor <- vector("list", length(relevant.Verticals.M))

  
  for(i in 1:length(relevant.Verticals.M)){ 
    # if(length(relevant.Verticals.M[[i]])>1){

   stream.Verticals.M <-  unlist(relevant.Verticals.M[[i]], recursive=F,use.names = F) 
   names(stream.Verticals.M) <- as.vector(unlist(lapply(relevant.Verticals.M[[i]],names)))
   stream.Verticals.G <-  unlist(relevant.Verticals.G, recursive=F,use.names = F) 
   names(stream.Verticals.G) <- as.vector(unlist(lapply(relevant.Verticals.G,names)))
    
    factors.M <- unique(unlist(stream.Verticals.M))
    factors.GT <- unique(unlist(stream.Verticals.G))
    differ <- setdiff(toupper(factors.GT) , toupper(factors.M))
    
    streamline.G <- lapply(stream.Verticals.G,function(y) y[which(!y%in% union(differ,tolower(differ)))]) 
    streamline.G <-lapply(streamline.G,function(y) y[length(y)!=1])
    streamline.G <- streamline.G[which(unlist(lapply(streamline.G,function(y) length(y)>1)))]
    streamline.G <- lapply(streamline.G, function(s) paste(s,collapse = "-"))
    # streamline.G <- split(streamline.G,names(streamline.G))
    
    # per.outcome <- vector("list", length(relevant.Verticals.M[[i]]))
    # for(k in 1:length(relevant.Verticals.M[[i]])){
    intersection <- vector("list", length(stream.Verticals.M))
    for(j in 1:length(stream.Verticals.M)){
      out.name <- names(stream.Verticals.M[j])
      rr <- paste(stream.Verticals.M[[j]], collapse = "-")
      qq <- expand.grid(list(rr), streamline.G[which(names(streamline.G)%in% out.name)], stringsAsFactors = F)
      intersection[[j]] <- any(apply(qq, 1, function(x) grepl(x[1],x[2])))
    }
    
    stream.Verticals.M <- stream.Verticals.M[which(unlist(intersection))]
    stream.Verticals.M <- split(stream.Verticals.M,names(stream.Verticals.M))
    stream.Verticals.M <- stream.Verticals.M[names(relevant.Verticals.M[[i]])]
    per.outcome <- vector("list", length(relevant.Verticals.M[[i]]))
    weight <- vector("list", length(relevant.Verticals.M[[i]]))
    for(j in 1:length(relevant.Verticals.M[[i]])){
      
      xx <- stream.Verticals.M[which(names(stream.Verticals.M)%in%names(relevant.Verticals.M[[i]][[j]]))]
      xx <- unlist(xx,recursive = F)
     per.outcome[[j]] <-  length(xx)/length(relevant.Verticals.M[[i]][[j]])
           weight[[j]]  <- length(relevant.Verticals.M[[i]][[j]])/length(unlist(relevant.Verticals.M[[i]], recursive = F))
    }
    # lapply(stream.Verticals.M, function(x), lenght(x)/length())
    

    
   
    
    # }
     # per.outcome.correctness <- as.list(mapply(function(x,y) sum(unlist(x))/length(y), x= per.outcome, y=relevant.Verticals.M[[i]]))
     names(per.outcome) <- attr(y$relevant.Sequences[[i]], "names")

     per.outcome<- as.data.frame(do.call(rbind, per.outcome)  )
     names(per.outcome) <- "score"
     per.outcome$weight <- unlist(weight) 
     weighted <-   weighted.mean(per.outcome[,1],per.outcome[,2],na.rm = T)
     per.outcome$weighted.mean <- weighted
     VertCor[[i]]  <- per.outcome
     
  # }else{VertCor[[i]] <-  NA }
    VertCor[[i]] <- list(VertCor[[i]])
   names(VertCor[[i]]) <- "Sequential.Correctness"
    attr(VertCor[[i]], "relevant.Sequences") <- relevant.Verticals.M[[i]]
    
    }


   
    # ======= Vertical Completeness


   VertCom <- vector("list", length(relevant.Verticals.M))
  for(i in 1:length(relevant.Verticals.M)){ 
    
   stream.Verticals.M <-  unlist(relevant.Verticals.M[[i]], recursive=F,use.names = F) 
   names(stream.Verticals.M) <- as.vector(unlist(lapply(relevant.Verticals.M[[i]],names)))
   stream.Verticals.G <-  unlist(relevant.Verticals.G, recursive=F,use.names = F) 
   names(stream.Verticals.G) <- as.vector(unlist(lapply(relevant.Verticals.G,names)))
    
    factors.M <- unique(unlist(stream.Verticals.M))
    factors.GT <- unique(unlist(stream.Verticals.G))
    differ <- setdiff(toupper(factors.GT) , toupper(factors.M))

    
    streamline.M <- lapply(stream.Verticals.M,function(y) y[which(!y%in% union(differ,tolower(differ)))]) 
    streamline.M <-lapply(streamline.M,function(y) y[length(y)!=1])
    streamline.M <- streamline.M[which(unlist(lapply(streamline.M,function(y) length(y)>1)))]
    streamline.M <- lapply(streamline.M, function(s) paste(s,collapse = "-"))
    # streamline.M <- split(streamline.M,names(streamline.M))
    
    # per.outcome <- vector("list", length(relevant.Verticals.G))
    # 
    # for(k in 1:length(relevant.Verticals.G)){
    # intersection <- vector("list", length(relevant.Verticals.G[[k]]))
    for(j in 1:length(stream.Verticals.G)){
      out.name <- names(stream.Verticals.G[j])
      rr <- paste(stream.Verticals.G[[j]], collapse = "-")
      qq <- expand.grid(list(rr), streamline.M[which(names(streamline.M)%in%out.name)], stringsAsFactors = F)
      inter <- apply(qq, 1, function(x) stri_extract_all_coll(x[1],x[2]))
      inter <- unlist(inter)
      if(!all(is.na(inter))){
      inter <- inter[!is.na(inter)]
      inter <- inter[which(nchar(inter)==max(nchar(inter)))][1]
      intersection[[j]] <- length(unlist(str_split(inter, "-")))/length(unlist(str_split(rr, "-")))
      }else{intersection[[j]] <- 0}
    }
    
    names(intersection) <- names(stream.Verticals.G)
    intersection <- split(intersection,names(intersection))
    weights <- lapply(intersection, function(x) length(x)/length(unlist(intersection)))

    # streamline.G <- split(streamline.G,names(streamline.G))
    per.outcome <- lapply(intersection, function(x) mean(unlist(x), na.rm = T))
      


     per.outcome<- as.data.frame(do.call(rbind, per.outcome)  )
     names(per.outcome) <- "score"
     per.outcome$weight <- unlist(weights) 
     weighted <-   weighted.mean(per.outcome[,1],per.outcome[,2],na.rm = T)
     per.outcome$weighted.mean <- weighted
     VertCom[[i]]  <- per.outcome
        
    
    VertCom[[i]] <- list(VertCom[[i]])
    

   names(VertCom[[i]]) <- "Sequential.Completeness"
    attr(VertCom[[i]], "relevant.Sequences") <- relevant.Verticals.G
    
  }
   
   
   
  # ====== FBeta Aggregation
  
    FBeta.Vert <- vector("list", length(relevant.Verticals.M))
    for(i in 1:length(FBeta.Vert)){
      FBeta <- (1+beta^2)*(as.numeric(VertCor[[i]]$Sequential.Correctness$weighted[1]) * as.numeric(VertCom[[i]]$Sequential.Completeness$weighted.mean[1]) )/
                  ((beta^2*as.numeric(VertCor[[i]]$Sequential.Correctness$weighted[1]))+ as.numeric(VertCom[[i]]$Sequential.Completeness$weighted.mean[1])) 
       if(is.na(FBeta)){FBeta <- 0}
      FBeta.Vert[[i]] <-  list(VertCor[[i]],VertCom[[i]],FBeta=FBeta)
      attr(FBeta.Vert[[i]], "model") <- attr(relevant.Verticals.M[[i]], "model4")
      attr(FBeta.Vert[[i]], "groundTruth") <- GT
      }  
   
###### Results 
if(any(class(x)=="character")){out <- vector("list", length(x))}
if(any(class(x)=="cna")){out <- vector("list", length(csfs.M))}
  for(i in 1:length(out)){
   
    if(agregate=="mean1"){

  # correctness
cor.sumLit <- unlist(attr(litCor[[i]], "relevant.Literals") ) |> length()
cor.sumCon <- unlist(str_split(unlist(attr(ConjCor[[i]], "relevant.Conjuncts") ), "[*]")) |> length()
cor.sumDis <- unlist(str_split(unlist(attr(DisjCor[[i]], "relevant.Disjunctions")), "[+,*]" )) |> length()
cor.sumVert <- unlist(attr(VertCor[[i]], "relevant.Sequences"), recursive = F)  |> length()
cor.total <- sum(cor.sumLit,cor.sumCon,cor.sumDis,cor.sumVert)
cor.weights <-  c(cor.sumLit/cor.total, cor.sumCon/cor.total, cor.sumDis/cor.total,cor.sumVert/cor.total)
cor.overall <- weighted.mean(c(litCor[[i]]$Literal.Correctness$weighted.mean[1], ConjCor[[i]]$Conjunction.Correctness$weighted.mean[1],
                      DisjCor[[i]]$Disjunction.Correctness$weighted.mean[1], as.numeric(VertCor[[i]]$Sequential.Correctness$weighted.mean[1])), cor.weights, na.rm = TRUE)
 # completeness
com.sumLit <- unlist(attr(litCom[[i]], "relevant.Literals") ) |> length()
com.sumCon <- unlist(str_split(unlist(attr(ConjCom[[i]], "relevant.Conjuncts") ), "[*]")) |> length()
com.sumDis <- unlist(str_split(unlist(attr(DisjCom[[i]], "relevant.Disjunctions") ), "[+,*]" )) |> length()
com.sumVert <- unlist(attr(VertCom[[i]], "relevant.Sequences"), recursive = F)  |> length()
com.total <- sum(com.sumLit,com.sumCon,com.sumDis,com.sumVert)
com.weights <-  c(com.sumLit/com.total, com.sumCon/com.total, com.sumDis/com.total,com.sumVert/com.total)
com.overall <- weighted.mean(c(litCom[[i]]$Literal.Completeness$weighted.mean[1], ConjCom[[i]]$Conjunction.Completeness$weighted.mean[1],
                      DisjCom[[i]]$Disjunction.Completeness$weighted.mean[1], as.numeric(VertCom[[i]]$Sequential.Completeness$weighted.mean[1])), com.weights, na.rm = TRUE)
  }

if(agregate=="mean2"){

  # correctness
cor.sumLit <- unlist(attr(litCor[[i]], "relevant.Literals") ) |> length()
cor.sumCon <- unlist(str_split(unlist(attr(ConjCor[[i]], "relevant.Conjuncts") ), "[*]")) |> length()
cor.sumDis <- unlist(str_split(unlist(attr(DisjCor[[i]], "relevant.Disjunctions")), "[+,*]" )) |> length()
cor.sumVert <- unlist(attr(VertCor[[i]], "relevant.Sequences"), recursive = F)  |> length()
# cor.total <- sum(cor.sumLit,cor.sumCon,cor.sumDis,cor.sumVert)
# cor.weights <-  c(cor.sumLit/cor.total, cor.sumCon/cor.total, cor.sumDis/cor.total,cor.sumVert/cor.total)
cor.overall <- mean(c(litCor[[i]]$Literal.Correctness$weighted.mean[1], ConjCor[[i]]$Conjunction.Correctness$weighted.mean[1],
                      DisjCor[[i]]$Disjunction.Correctness$weighted.mean[1], as.numeric(VertCor[[i]]$Sequential.Correctness$weighted.mean[1])), na.rm = TRUE)
 # completeness
com.sumLit <- unlist(attr(litCom[[i]], "relevant.Literals") ) |> length()
com.sumCon <- unlist(str_split(unlist(attr(ConjCom[[i]], "relevant.Conjuncts") ), "[*]")) |> length()
com.sumDis <- unlist(str_split(unlist(attr(DisjCom[[i]], "relevant.Disjunctions") ), "[+,*]" )) |> length()
com.sumVert <- unlist(attr(VertCom[[i]], "relevant.Sequences"), recursive = F)  |> length()
# com.total <- sum(com.sumLit,com.sumCon,com.sumDis,com.sumVert)
# com.weights <-  c(com.sumLit/com.total, com.sumCon/com.total, com.sumDis/com.total,com.sumVert/com.total)
com.overall <- mean(c(litCom[[i]]$Literal.Completeness$weighted.mean[1], ConjCom[[i]]$Conjunction.Completeness$weighted.mean[1],
                      DisjCom[[i]]$Disjunction.Completeness$weighted.mean[1], as.numeric(VertCom[[i]]$Sequential.Completeness$weighted.mean[1])), na.rm = TRUE)
}
    

if(agregate=="mean2"){

  # correctness
cor.sumLit <- unlist(attr(litCor[[i]], "relevant.Literals") ) |> length()
cor.sumCon <- unlist(str_split(unlist(attr(ConjCor[[i]], "relevant.Conjuncts") ), "[*]")) |> length()
cor.sumDis <- unlist(str_split(unlist(attr(DisjCor[[i]], "relevant.Disjunctions")), "[+,*]" )) |> length()
cor.sumVert <- unlist(attr(VertCor[[i]], "relevant.Sequences"), recursive = F)  |> length()
# cor.total <- sum(cor.sumLit,cor.sumCon,cor.sumDis,cor.sumVert)
# cor.weights <-  c(cor.sumLit/cor.total, cor.sumCon/cor.total, cor.sumDis/cor.total,cor.sumVert/cor.total)
cor.overall <- mean(c(litCor[[i]]$Literal.Correctness$weighted.mean[1], ConjCor[[i]]$Conjunction.Correctness$weighted.mean[1],
                      DisjCor[[i]]$Disjunction.Correctness$weighted.mean[1], as.numeric(VertCor[[i]]$Sequential.Correctness$weighted.mean[1])), na.rm = TRUE)
 # completeness
com.sumLit <- unlist(attr(litCom[[i]], "relevant.Literals") ) |> length()
com.sumCon <- unlist(str_split(unlist(attr(ConjCom[[i]], "relevant.Conjuncts") ), "[*]")) |> length()
com.sumDis <- unlist(str_split(unlist(attr(DisjCom[[i]], "relevant.Disjunctions") ), "[+,*]" )) |> length()
com.sumVert <- unlist(attr(VertCom[[i]], "relevant.Sequences"), recursive = F)  |> length()
# com.total <- sum(com.sumLit,com.sumCon,com.sumDis,com.sumVert)
# com.weights <-  c(com.sumLit/com.total, com.sumCon/com.total, com.sumDis/com.total,com.sumVert/com.total)
com.overall <- mean(c(litCom[[i]]$Literal.Completeness$weighted.mean[1], ConjCom[[i]]$Conjunction.Completeness$weighted.mean[1],
                      DisjCom[[i]]$Disjunction.Completeness$weighted.mean[1], as.numeric(VertCom[[i]]$Sequential.Completeness$weighted.mean[1])), na.rm = TRUE)
  }    


  if(agregate=="harmonic1"){

  # correctness
cor.overall <- harmonic.mean(c(litCor[[i]]$Literal.Correctness$weighted.mean[1], ConjCor[[i]]$Conjunction.Correctness$weighted.mean[1],
                      DisjCor[[i]]$Disjunction.Correctness$weighted.mean[1], as.numeric(VertCor[[i]]$Sequential.Correctness$weighted.mean[1])), na.rm = TRUE, zero=FALSE)
if(is.na(cor.overall))cor.overall <- 0
 # completeness
com.overall <- harmonic.mean(c(litCom[[i]]$Literal.Completeness$weighted.mean[1], ConjCom[[i]]$Conjunction.Completeness$weighted.mean[1],
                      DisjCom[[i]]$Disjunction.Completeness$weighted.mean[1], as.numeric(VertCom[[i]]$Sequential.Completeness$weighted.mean[1])), na.rm = TRUE,zero=FALSE)
if(is.na(com.overall))com.overall <- 0
  }
    
    
    
    if(agregate=="harmonic2"){

  # correctness
cor.sumLit <- unlist(attr(litCor[[i]], "relevant.Literals") ) |> length()
cor.sumCon <- unlist(str_split(unlist(attr(ConjCor[[i]], "relevant.Conjuncts") ), "[*]")) |> length()
cor.sumDis <- unlist(str_split(unlist(attr(DisjCor[[i]], "relevant.Disjunctions")), "[+,*]" )) |> length()
cor.sumVert <- unlist(attr(VertCor[[i]], "relevant.Sequences"), recursive = F)  |> length()
cor.total <- sum(cor.sumLit,cor.sumCon,cor.sumDis,cor.sumVert)
cor.weights <-  c(cor.sumLit/cor.total, cor.sumCon/cor.total, cor.sumDis/cor.total,cor.sumVert/cor.total)
cor.values <- c(litCor[[i]]$Literal.Correctness$weighted.mean[1], ConjCor[[i]]$Conjunction.Correctness$weighted.mean[1],
                      DisjCor[[i]]$Disjunction.Correctness$weighted.mean[1], as.numeric(VertCor[[i]]$Sequential.Correctness$weighted.mean[1]))
cor.overall <- sum(cor.weights) / sum(cor.weights / cor.values)
 # completeness
com.sumLit <- unlist(attr(litCom[[i]], "relevant.Literals") ) |> length()
com.sumCon <- unlist(str_split(unlist(attr(ConjCom[[i]], "relevant.Conjuncts") ), "[*]")) |> length()
com.sumDis <- unlist(str_split(unlist(attr(DisjCom[[i]], "relevant.Disjunctions") ), "[+,*]" )) |> length()
com.sumVert <- unlist(attr(VertCom[[i]], "relevant.Sequences"), recursive = F)  |> length()
com.total <- sum(com.sumLit,com.sumCon,com.sumDis,com.sumVert)
com.weights <-  c(com.sumLit/com.total, com.sumCon/com.total, com.sumDis/com.total,com.sumVert/com.total)
com.values <- c(litCom[[i]]$Literal.Completeness$weighted.mean[1], ConjCom[[i]]$Conjunction.Completeness$weighted.mean[1],
                      DisjCom[[i]]$Disjunction.Completeness$weighted.mean[1], as.numeric(VertCom[[i]]$Sequential.Completeness$weighted.mean[1]))
com.overall <- sum(com.weights) / sum(com.weights / com.values)
  }    
    

FBeta.overall  <- (1+beta^2)*(cor.overall * com.overall )/((beta^2*cor.overall)+ com.overall)
if(is.na(FBeta.overall)){FBeta.overall <- 0}
  

    scoreboard <-  matrix(nrow = 3,ncol=5, dimnames = list(c("Correctness", "Completeness","FBeta"),
                                        c("Literals", "Conjunctions", "Disjunctions", "Sequences","Weighted.Total")))
    scoreboard[1,] <- c(litCor[[i]]$Literal.Correctness$weighted.mean[1],
               ConjCor[[i]]$Conjunction.Correctness$weighted.mean[1],
               DisjCor[[i]]$Disjunction.Correctness$weighted.mean[1],
               as.numeric(VertCor[[i]]$Sequential.Correctness$weighted.mean[1]), cor.overall) 
     scoreboard[2,] <- c(litCom[[i]]$Literal.Completeness$weighted.mean[1],
               ConjCom[[i]]$Conjunction.Completeness$weighted.mean[1],
               DisjCom[[i]]$Disjunction.Completeness$weighted.mean[1],
               as.numeric(VertCom[[i]]$Sequential.Completeness$weighted.mean[1]), com.overall) 
     scoreboard[3,] <- c(FBeta.Lit[[i]]$FBeta,
               FBeta.Conj[[i]]$FBeta,
               FBeta.Disj[[i]]$FBeta,
               as.numeric(FBeta.Vert[[i]]$FBeta), FBeta.overall) 


out[[i]] <- list(correctness=list(litCor[[i]], ConjCor[[i]],DisjCor[[i]],VertCor[[i]]),
                 completeness=list(litCom[[i]], ConjCom[[i]],DisjCom[[i]],VertCom[[i]]),
                 FBeta=list(FBeta.Lit[[i]], FBeta.Conj[[i]],FBeta.Disj[[i]],FBeta.Vert[[i]]), scoreboard=scoreboard) 

}  


out <- structure(out, class = "quality")
out
}


print.quality <- function(x, ...) {
  out <- vector("list",length(x))
  for(i in 1:length(x)){
    # xx <- unlist(x[[i]],recursive = F)

    out[[i]] <-  list(Ground.Truth=attr(x[[i]]$FBeta[[1]], "groundTruth"), 
                      model=attr(x[[i]]$FBeta[[1]], "model"), x[[i]][[4]] )
  }
  print(out)
}

fast_isSubmodel <- function(x, y){
  a1 <- hstrsplit(lhs(x), c("+", "*"))
  a2 <- hstrsplit(lhs(y), c("+", "*"))
  nms <- unique(c(unlist(a1), unlist(a2)))
  b1 <- rapply(a1, match, table = nms, how = "replace")
  b2 <- rapply(a2, match, table = nms, how = "replace")
  cna:::C_is_submodel(b1, b2[[1]], strict = FALSE)
}  

asf.overlap4 <- function (x, score.only = FALSE, unique = T) 
{
  
  x <- noblanks(x)
  if (length(x) == 1&&score.only==F&&unique==F) {
    return(list(max.intersect=lhs(x),common.leaves=unlist(str_extract_all(lhs(x), "[A-Za-z]")),
                score=length(unlist(str_extract_all(lhs(x), "[*,+]"))) + length(unlist(str_extract_all(lhs(x), "[A-Za-z]")))))
  }
  if (length(x) == 1&&score.only==F&&unique==T) {
    return(list(max.intersect=lhs(x),common.leaves=unique(unlist(str_extract_all(lhs(x), "[A-Za-z]"))),
                score=length(unlist(str_extract_all(lhs(x), "[*,+]"))) + length(unique(unlist(str_extract_all(lhs(x), "[A-Za-z]"))))))
  }
  if (length(x) == 1&&score.only==T&&unique==T) {
    return(length(unlist(str_extract_all(lhs(x), "[*,+]"))) + length(unique(unlist(str_extract_all(lhs(x), "[A-Za-z]")))))
  }
  if (length(x) == 1&&score.only==T&&unique==F) {
    return(length(unlist(str_extract_all(lhs(x), "[*,+]"))) + length(unlist(str_extract_all(lhs(x), "[A-Za-z]"))))
  }
  if (!all(condType(x) %in% "stdAtomic")) {
    stop("x must be a character vector of asf only")
  }
  out <-   unique(rhs(x)) 
  if (length(out)>1) {
    stop("All asf in x must have the same outcome")
  }
  
  
  y <-  as.data.frame(cbind(unlist(x),getComplexity(unlist(x))))
  colnames(y) <- c("models", "complexity")     
  y$complexity <- as.numeric(y$complexity)
  y <- y[order(y$complexity, decreasing = F),]
  # y <- split(y,y$complexity)
  
  # identify ground (starting point for submodel search, models with smallest complexity)
  # dynamic determination of the amount of models used to which asf.overlap.brute is applied
  max.comp <- y[1:min(nrow(y),4),]
  max.comp <- ifelse(max(max.comp$complexity)<5, 4, ifelse(max(max.comp$complexity)<9, 3, 2))
  ground <- y[1:min(nrow(y),max.comp),]$models
  # cap ground at 5
  # if(length(ground)>2)ground <- ground[1:2]
  # identify max.intersect of ground (with brute force)
  upper.bound <- asf.overlap.brute(as.list(ground))
  # stop if ground is everything in x
  rest <- setdiff(x,ground)
  if((length(rest)==0||upper.bound$score==0)&&score.only==F) return(upper.bound)
  if((length(rest)==0||upper.bound$score==0)&&score.only==T) return(upper.bound$score)
  # if ground is not everything in x, prepare submodel test with rest
  stop <- FALSE
  test <- as.list(paste0(upper.bound$max.intersect, "<->", out))
  while(stop==FALSE){
    
    result <- lapply(test, function(x) lapply(rest, function(y) is.submodel(x, y)))
    r <- lapply(result, function(x) all(unlist(x))) 
    if(any(unlist(r))){stop <- TRUE
    }else{
      
      # Disect test
      f <- lhs(unlist(test)) 
      dis <- str_split(f,"\\+")
      conj.test <- lapply(dis, function(y) str_split(y, "\\*"))
      
      test <- vector("list", length(conj.test))
      for(i in 1:length(conj.test)){
        xx <- vector("list", length(conj.test[[i]]))
        for(j in 1:length(conj.test[[i]])){
          xxx <- vector("list", length(conj.test[[i]][[j]]))
          for(k in 1:length(conj.test[[i]][[j]])){
            ma <- conj.test[[i]]
            ma[[j]] <- ma[[j]][-k]
            ma <- lapply(ma, function(x) paste(x,collapse="*"))
            ma <- lapply(ma, function(x)if(x==""){x <- NULL}else{x})
            if(length(unlist(ma)>0)){ma <- noblanks(paste(paste(unlist(ma),collapse = "+", recycle0=T),"<->",out))
            }else{ma <- NULL}
            xxx[[k]] <- ma
          }
          xx[[j]] <- xxx   
        }
        test[[i]] <- xx  
      }
      test <- as.list(unique(unlist(test)))  
      if(length(test)==0)stop <- TRUE
    }
  } 
  # one max.intersect, among possibly many: 
  max.intersect <- unique(test[unlist(r)]) 
  # find common leaves
  all.dnf <- lhs(unlist(x))
  all.leaves <- str_extract_all(all.dnf, "[A-Za-z]")
  if(unique==T){common.leaves <- unique(Reduce(intersect, all.leaves ))
  }else{common.leaves <- Reduce(intersect, all.leaves )}
  # calculate overlap score
  score.inter <- length(common.leaves)
  score.connect <- if(length(max.intersect)>0){length(unlist(str_extract_all(max.intersect[1], "[*,+]")))}else{0}
  score <- score.inter+score.connect
  # return results    
  if(score.only==FALSE){
    return(list(max.intersect=unlist(max.intersect), common.leaves=common.leaves, score=score))
  }else{return(score)  
    
  }
}

asf.overlap.brute <- function (x, all = FALSE, score.only = FALSE) {
     if (length(x) == 1) {
        return(list(max.intersect=lhs(x),common.leaves=str_extract_all(lhs(x), "[a-zA-Z]"),
                    score=length(unlist(str_extract_all(lhs(x), "[*,+]"))) + length(unlist(str_extract_all(lhs(x), "[A-Za-z]")))))
     }
     # isolate lhs and standardize syntax
    xx <- stdCond(lhs(x))
    dis <- str_split(xx,"\\+")
    max.overlap.dis <- min(unlist(lapply(dis,length)))
    conj <- lapply(dis, function(y) str_split(y, "\\*"))
    xx.grid <- expand.grid(conj,stringsAsFactors = F)
    xx.grid <- as.data.frame(xx.grid)
    attr(xx.grid, "out.attrs") <- NULL
    # str(xx.grid[1,])
    maxlen <- max(sapply(xx.grid,length))
    xx.inter <- lapply(seq(maxlen),function(i) Reduce(intersect,lapply(xx.grid,"[[",i)))
    # lapply(xx.inter, function(x) if(identical(x, character(0)))0 else x)
    if(length(unlist(xx.inter))==0&&score.only==FALSE)return(list(max.intersect=NULL, common.leaves=NULL, score=0))
    if(length(unlist(xx.inter))==0&&score.only==TRUE)return(0)
    xx.grid$overlap <- xx.inter
    xx.grid$length <- apply(xx.grid,1, function(y) length(y$overlap))
    xx.select <- xx.grid[which(xx.grid$length>0),]
    xx.select <- xx.select[order(xx.select$length),]
      # build overlaps
    b <- xx.select$overlap 
    b <- lapply(b, function(x) paste(x, collapse="*"))
    b <- unlist(b)

    if(length(b)>1){
    tt <- xx.select[,-which(names(xx.select)%in%c("overlap","length"))]
    n <- min(length(b), max.overlap.dis)
    xx.overlap <- vector("list",n)
    if(score.only==FALSE){
      for(i in 1:n){
    xx.overlap[[i]] <- combn(1:length(xx.select$overlap), i ,simplify=F)
      }
    }else{
      for(i in n){
    xx.overlap[[i]] <- combn(1:length(xx.select$overlap), i ,simplify=F)
      }
    }
    xx.overlap <- unlist(xx.overlap,recursive = F)
    xx.overlap <- xx.overlap[unlist(lapply(xx.overlap, function(x)!any(apply(tt[x,], 2, function(x)duplicated(x)))))]
    
    if(length(xx.overlap)>0){
      overlaps <- lapply(xx.overlap,function(x) paste(b[unlist(x)], collapse="+"))
      overlaps <- lapply(overlaps, stdCond)
      com <- max(nchar(unlist(overlaps)))
      overlaps <- c(b[nchar(b)==com],overlaps)
      }else{overlaps <- b}
    overlaps <- unique(unlist(overlaps))
    }else{overlaps <- b}
    max.overlaps <- overlaps[which(nchar(overlaps)==max(nchar(overlaps)))] 
    
    # Calculate score
    cl <- lapply(xx, function(y) str_extract_all(y, "[a-zA-Z]"))
    cl <- unlist(cl,recursive = F)
    common.leaves <- Reduce(intersect,cl)

    # common.leaves <- str_extract_all(max.overlaps, "[a-zA-Z]")
    
    combinations <- vector("list", length(common.leaves))
    for(j in 1:length(common.leaves)){
      combinations[[j]] <- combn(1:length(common.leaves), j ,simplify=F)
      
    }
    
    # base <- lapply(x[[i]], function(s) asf.overlap3(s,score.only = T))
    base <- length(unlist(common.leaves))
    # base <- sum(unlist(base))
    if(length(combinations)>1){
      r <- split(2:length(combinations), 2:length(combinations) %% 2 == 0)
      # add odd overlaps
      odd <- r[["FALSE"]]
      toAdd <- lapply(combinations[odd], function(y) lapply(y, function(s) common.leaves[[1]][common.leaves[[1]]%in% Reduce(intersect, common.leaves[s] )]))
      toAdd <- length(unlist(toAdd))
      
      # subtract even overlaps
      even <- r[["TRUE"]]
      toSubtract <- lapply(combinations[even], function(y) lapply(y, function(s) common.leaves[[1]][common.leaves[[1]]%in% Reduce(intersect, common.leaves[s] )]))
      toSubtract <-length(unlist(toSubtract))      
      score.inter <- base+toAdd-toSubtract
    }else{score.inter <- base}

    max.intersect <- unique(stdCond2(max.overlaps))
    # score.inter <- length(common.leaves)
    score.connect <- length(unlist(str_extract_all(max.intersect[1], "[*,+]")))
    score <- score.inter+score.connect
    # output
    if(score.only==FALSE){
    if(all==FALSE){return(list(max.intersect=max.intersect , common.leaves=common.leaves, score=score))}
    if(all==TRUE){return(list(all.intersect=unique(stdCond2(overlaps)), common.leaves=common.leaves, score=score))}
    }else{
     return(score) 
    }
}

condType <- function(cond, x = full.ct(cond)){
  out <- cna:::getCondType(cond, cna:::ctInfo(configTable(x)))
  as.vector(out)
}


stdCond2 <- function (x) 
{
    l <- strsplit(noblanks(x), "+", fixed = TRUE)
    if (length(l) == 0L) 
        return(character(0))
    # l <- lapply(l, unique.default)
    u <- unlist(l, use.names = FALSE, recursive = FALSE)
    out <- cna:::C_relist_Char(cna:::C_mconcat(strsplit(u, "*", fixed = TRUE), 
        "*", sorted = TRUE), lengths(l))
    cna:::C_mconcat(out, "+", sorted = TRUE)
}

switch.case <- function(x){
ifelse(useful::find.case(x), tolower(x),toupper(x))
}

negate3 <- function(x){
  x <- str_split(x, "[+]")
  x <- lapply(x, function(r) str_split(r, "[*]")) 
  x <- lapply(x, function(s) lapply(s, function(r) switch.case(r))) 
  # x <- lapply(x, function(y)paste0("!(",y,")"))
   # x <- lapply(x, function(y) lapply(y,function(d)rreduce(getCond(selectCases(d)))))
  x <- lapply(x, function(s) lapply(s, function(r) paste(r, collapse = "+"))) 
  x <- lapply(x,function(x) lapply(x, stdCond))
  x <- lapply(x, function(y) lapply(y, function(r) paste0("(",r,")")))
  x <- unlist(x)
   # x <- paste0("(",x,")")
   x <- paste(x,collapse="*")
   out <- factor.out2(x)
   noblanks(out)
}

factor.out2 <- function(x){
dis <- str_split(x, "\\++(?=[^()]*(?:\\(|$))")
conj <- lapply(dis, function(x) str_split(x, "\\*+(?=[^()]*(?:\\(|$))"))
conj <- lapply(conj, function(x)lapply(x, function(y) str_remove_all(y,"\\(")))
conj <- lapply(conj, function(x)lapply(x, function(y) str_remove_all(y,"\\)")))
dis_in <- lapply(conj, function(x)lapply(x, function(y) str_split(y,"\\+")))
dis_in <- lapply(dis_in, function(x) lapply(x,function(y) expand.grid(y)))
# apply(dis_in[[1]][[1]],1,
dis_red <- lapply(dis_in, function(x) lapply(x, function(q) apply(q, 1, function(e) str_split(e, "\\*", simplify = F))))
dis_red <- lapply(dis_red, function(x) lapply(x, function(q) lapply(q, unlist)))
# unlist(dis_red, recursive = F) |> unlist(recursive = F)
dis_red <- lapply(dis_red,function(s) lapply(s, function(w) lapply(w, function(r) unique(r))))

dis_red <- lapply(dis_red, function(w) lapply(w, function(x) lapply(x, function(s) if(any(duplicated(toupper(s)))){s <- NULL}else{s})))
dis_red <- lapply(dis_red, function(x) lapply(x, function(w) Filter(Negate(is.null), w)))
dis_red <- lapply(dis_red, function(w) lapply(w, function(q) lapply(q, function(x) paste(x, collapse = "*"))))
dis_red <- lapply(dis_red, function(s) lapply(s, function(r) lapply(r, cna:::stdCond)))
dis_red <- unlist(dis_red)
conj <- dis_red[!duplicated(dis_red)]
out <- paste(noblanks(conj), collapse = "+")
out
}


LCR <- function(x,GT){
  out <- is.submodel(x,GT)
  as.numeric(out)
  }

get_submodels <- function(solution) {
    outcome <- rhs(solution)
    out <- QCApro::submodels(solution)$submodels
    out <- out[out != '']
    out <- paste0(out,"<->",outcome) |> noblanks()
    return(out)
}


wrongness <- function(M,GT) {
     x <- get_submodels(M)
     y <- get_submodels(GT)
     sum(!x %in% y) / length(x)
     
}

NCR <- function(M,GT) {
     x <- get_submodels(M)
     y <- get_submodels(GT)
     sum(x %in% y) / length(x)
     
}





