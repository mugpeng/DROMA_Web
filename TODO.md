# TODO

## Major

- [x] Add PDO drug and rna
- [x] Reonline
- [x] change to z-score 
- [x] server check1 (process die)
- [ ] single data compare, find cell line, compare types(CXP)
  - [ ] Try Co_Occurence(like maftools) part for drugs
- [ ] Add enrichment methods
- [ ] Harmonized drug and cell names
  - [ ] how to build a name checker agent?

- [ ] add raw dose viability part, recomputeAUC using PharmacoGx method,
- [ ] How to allow user analysis data without loading all datasets locally
- [ ] reload normalization will cause process dead
- [ ] Add combintation drugs results parts
  - [ ] in house
  - [ ] synergy prediction
- [ ] server check2 (make shiny load faster)
  - Preloading noticed screen
  - options(shiny.idle_timeout = 0)
  - Using the future and promises packages to run computations asynchronously
- [ ] Add data ranking check(allow user highlight interested features), drug rank(select a cell) and cell rank(select a feature)
- [ ] Add PDO WES
- [ ] add chemical structure info 
- [ ] https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-014-0367-1/figures/1 seems better to use rank based method (how to solve the direction problem, seems ok, smaller means another direction, and the median rank is to devide into two parts)
  - try to see the performance of rank method
- [ ] more detailed cellline data
  - [ ] [Cellosaurus cell line MDA-MB-361 (CVCL_0620)](https://www.cellosaurus.org/CVCL_0620)




## Minor

- [x] filter continous data(zero SD)
- [x] add choice to download csv data
- [ ] what para orchestra official used for calculate arc
- [ ] recheck cellline WES results
- [ ] Add drug annotation for drug screen in batch mode
- [ ] volcano plot seems a bug wrongly highlight top5
- [ ] Add compare methods
- [ ] Add user counts
- [ ] adapt to mobile user



## Function

- [ ] Add ecDNA predict method
  - MYC ecDNA promotes intratumour heterogeneity and plasticity in PDAC, AmpliconArchitect
- [ ] 



# Similar products

TDC:

![image-20250328114629463](images/image-20250328114629463.png)

datasets for Txgemma/txagent, Marinka Zitnik



# datasets

发出来了吗？

![image-20250322182104830](images/image-20250322182104830.png)



