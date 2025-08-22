[![License: GPL](https://img.shields.io/badge/License-GPL-yellow)](https://github.com/AspirinCode/awesome-AI4MolConformation-MD)

## awesome-AI4MolConformation-MD
List of **molecules ( small molecules, RNA, peptide, protein, enzymes, antibody, and PPIs) conformations** and **molecular dynamics (force fields)** using **generative artificial intelligence** and **deep learning**


![Protein Space and Conformations](https://github.com/AspirinCode/awesome-AI4MolConformation-MD/blob/main/figure/afpro.png)


**Updating ...**  


## Menu

  - [Deep Learning-molecular conformations](#deep-learning-molecular-conformations)

| Menu | Menu | Menu | Menu |
| ------ | :---------- | ------ | ------ |
| [Reviews](#reviews) | [Datasets and Package](#datasets-and-package) | [Molecular dynamics](#molecular-dynamics) | [Molecular Force Fields](#molecular-force-fields) |
| [MD Engines-Frameworks](#md-engines-frameworks) | [AI4MD Engines-Frameworks](#ai4md-engines-frameworks) | [MD Trajectory Processing-Analysis](#md-trajectory-processing-analysis) | [CGMD](#cgmd) |
| [AI4MD](#ai4md) | [Neural Network Potentials](#neural-network-potentials) | [Free Energy Perturbation](#free-energy-perturbation) | [Ab Initio](#ab-initio) |
| [Neural Molecular Force Fields](#neural-molecular-force-fields) | [Neural Reactive Potential](#neural-reactive-potential) | [Reactive Force Fields](#reactive-force-fields) |  |
| [AlphaFold-based](#alphaFold-based) | [Autoregressive-Based](#autoregressive-Based)  | [LSTM-based](#lstm-based) | [Transformer-based](#transformer-based) |
| [VAE-based](#vae-based) | [GAN-based](#gan-based) | [Flow-based](#flow-based) | [Diffusion-based](#diffusion-based) |
| [Score-Based](#score-Based) | [Energy-based](#energy-based) | [Bayesian-based](#bayesian-based) | [Active Learning-based](#active-learning-based) |
| [GNN-based](#gnn-based) | [LLM-MD](#llm-md) |  |  |



  - [Molecular conformational ensembles by methods](#molecular-conformational-ensembles-by-methods)

| Menu | Menu | Menu |
| ------ | :---------- | ------ |
| [Small molecule conformational ensembles](#small-molecule-conformational-ensembles) | [RNA conformational ensembles](#rna-conformational-ensembles) | [Peptide conformational ensembles](#peptide-conformational-ensembles) |
| [Protein conformational ensembles](#protein-conformational-ensembles) | [Enzymes conformational ensembles](#enzymes-conformational-ensembles)  | [Antibody conformational ensembles](#antibody-conformational-ensembles) |
| [Ligand-Protein conformational ensembles](#ligand-protein-conformational-ensembles) | [RNA-Peptide conformational ensembles](#rna-peptide-conformational-ensembles) | [PPI conformational ensembles](#ppi-conformational-ensembles) |
| [Antibody-Protein conformational ensembles](#antibody-protein-conformational-ensembles) | [Nucleic acid-Protein conformational ensembles](#nucleic-acid-protein-conformational-ensembles) | [Material ensembles](#material-ensembles) |
| [Nucleic acid-Ligand conformational ensembles](#nucleic-acid-ligand-conformational-ensembles) |  |  |





## Reviews



* **A critical review of machine learning interatomic potentials and Hamiltonian** [2025]  
Li, Y.; Zhang, X.; Liu, M.; Shen, L.  
  [J. Mater. Inf. (2025)](http://dx.doi.org/10.20517/jmi.2025.17)  

* **Generation of protein dynamics by machine learning** [2025]  
Janson, Giacomo, and Michael Feig.  
  [Current Opinion in Structural Biology 93 (2025)](https://doi.org/10.1016/j.sbi.2025.103115)  

* **Beyond static structures: protein dynamic conformations modeling in the post-AlphaFold era** [2025]  
 Xinyue Cui, Lingyu Ge, Xia Chen, Zexin Lv, Suhui Wang, Xiaogen Zhou, Guijun Zhang.  
  [Briefings in Bioinformatics (2025)](https://doi.org/10.1093/bib/bbaf340)  

* **From sequence to protein structure and conformational dynamics with artificial intelligence/machine learning** [2025]  
 Alexander M. Ille, Emily Anas, Michael B. Mathews, Stephen K. Burley.  
  [Struct. Dyn. 12, 030902 (2025)](https://doi.org/10.1063/4.0000765)  

* **Application of machine learning interatomic potentials in heterogeneous catalysis** [2025]  
 Olajide, Gbolagade, Khagendra Baral, Sophia Ezendu, Ademola Soyemi, and Tibor Szilvasi.  
  [Journal of Catalysis (2025)](https://doi.org/10.1016/j.jcat.2025.116202)  

* **The evolution of machine learning potentials for molecules, reactions and materials** [2025]   
 Xia, Junfan and Zhang, Yaolong and Jiang, Bin.  
  [Chem. Soc. Rev. (2025)](https://doi.org/10.1039/D5CS00104H)  

* **Advancing Molecular Simulations: Merging Physical Models, Experiments, and AI to Tackle Multiscale Complexity** [2025]   
 Giorgio Bonollo, Gauthier Trèves, Denis Komarov, Samman Mansoor, Elisabetta Moroni, and Giorgio Colombo.  
  [J. Phys. Chem. Lett. (2025)](https://doi.org/10.1021/acs.jpclett.5c00652)  

* **A comparison of probabilistic generative frameworks for molecular simulations** [2025]   
 Richard John, Lukas Herron, Pratyush Tiwary.  
  [J. Chem. Phys. (2025)](https://doi.org/10.1063/5.0249683)  

* **Recent Advances in Machine Learning and Coarse-Grained Potentials for Biomolecular Simulations and Their Applications** [2025]   
 B. Poma A, Hinostroza Caldas A, Cofas-Vargas L, Jones M, L. Ferguson A, Medrano Sandonas L.  
  [JChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-vxjlk)  

* **Recent Advances in Simulation Software and Force Fields: Their Importance in Theoretical and Computational Chemistry and Biophysics** [2024]   
 Christophe Chipot.  
  [J. Phys. Chem. B (2024)](https://doi.org/10.1021/acs.jpcb.4c06231)  

* **Graph theory approaches for molecular dynamics simulations** [2024]   
 Patel AC, Sinha S, Palermo G.  
  [Quarterly Reviews of Biophysics. (2024)](https://doi.org/10.1017/S0033583524000143)  

* **Deep learning for intrinsically disordered proteins: From improved predictions to deciphering conformational ensembles** [2024]   
 Erdős, G., & Dosztányi, Z.   
  [Current opinion in structural biology (2024)](https://doi.org/10.1016/j.sbi.2024.102950)  

* **Recent advances in protein conformation sampling by combining machine learning with molecular simulation** [2024]   
 Tang, Y., Yang, Z., Yao, Y., Zhou, Y., Tan, Y., Wang, Z., Pan, T., Xiong, R., Sun, J. and Wei, G.   
  [Chinese Physics B. (2024)](https://iopscience.iop.org/article/10.1088/1674-1056/ad1a92)  

* **Alchemical Transformations and Beyond: Recent Advances and Real-World Applications of Free Energy Calculations in Drug Discovery** [2024]   
 Qian, Runtong, Jing Xue, You Xu, and Jing Huang.   
  [J. Chem. Inf. Model. (2024)](https://doi.org/10.1021/acs.jcim.4c01024)  

* **The need to implement FAIR principles in biomolecular simulations** [2024]   
 Amaro, Rommie, Johan Åqvist, Ivet Bahar, Federica Battistini, Adam Bellaiche, Daniel Beltran, Philip C. Biggin et al.   
  [arXiv:2407.16584 (2024)](https://arxiv.org/abs/2407.16584)  

* **An overview about neural networks potentials in molecular dynamics simulation** [2024]   
 Martin‐Barrios, Raidel, Edisel Navas‐Conyedo, Xuyi Zhang, Yunwei Chen, and Jorge Gulín‐González.   
  [International Journal of Quantum Chemistry 124.11 (2024)](https://doi.org/10.1002/qua.27389)  

* **Artificial Intelligence Enhanced Molecular Simulations** [2023]   
 Zhang, Jun, Dechin Chen, Yijie Xia, Yu-Peng Huang, Xiaohan Lin, Xu Han, Ningxi Ni et al.   
  [J. Chem. Theory Comput. (2023)](https://doi.org/10.1021/acs.jctc.3c00214)  

* **Machine Learning Generation of Dynamic Protein Conformational Ensembles** [2023]   
 Zheng, Li-E., Shrishti Barethiya, Erik Nordquist, and Jianhan Chen.   
  [Molecules 28.10 (2023)](https://doi.org/10.3390/molecules28104047)  





## Datasets and Package


### Datasets




* **DynaRepo: The repository of macromolecular conformational dynamics** [2025]  
Omid Mokhtari, Emmanuelle Bignon, Hamed Khakzad, Yasaman Karami.   
[bioRxiv (2025)](https://doi.org/10.1101/2025.08.14.670260) | [data](https://dynarepo.inria.fr/)  

* **MS25: Materials Science-Focused Benchmark Data Set for Machine Learning Interatomic Potentials** [2025]  
Tristan Maxson, Ademola Soyemi, Xinglong Zhang, Benjamin W. J. Chen, and Tibor Szilvási.   
[J. Chem. Inf. Model. (2025)](https://doi.org/10.1021/acs.jcim.5c01262) | [data](https://doi.org/10.5281/zenodo/10901820)  

* **A Benchmark for Quantum Chemistry Relaxations via Machine Learning Interatomic Potentials** [2025]  
Fu, Cong, Yuchao Lin, Zachary Krueger, Wendi Yu, Xiaoning Qian, Byung-Jun Yoon, Raymundo Arróyave et al.  
[arXiv:2506.23008 (2025)](https://doi.org/10.48550/arXiv.2506.23008) | [data](https://huggingface.co/divelab)  

* **The Open Molecules 2025 (OMol25) Dataset, Evaluations, and Models** [2025]  
Levine, Daniel S., Muhammed Shuaibi, Evan Walter Clark Spotte-Smith, Michael G. Taylor, Muhammad R. Hasyim, Kyle Michel, Ilyes Batatia, G'abor Cs'anyi, Misko Dzamba, Peter K. Eastman, Nathan C. Frey, Xiang Fu, Vahe Gharakhanyan, Aditi S. Krishnapriyan, Joshua A. Rackers, Sanjeev Raja, Ammar Rizvi, Andrew S. Rosen, Zachary W. Ulissi, Santiago Vargas, C. Lawrence Zitnick, Samuel M. Blau and Brandon M. Wood.   
[arXiv:2505.08762 (2025)](https://arxiv.org/abs/2505.08762) | [data](https://huggingface.co/facebook/OMol25)  

* **Enhanced Sampling, Public Dataset and Generative Model for Drug-Protein Dissociation Dynamics** [2025]  
Maodong Li, Jiying Zhang, Bin Feng, Wenqi Zeng, Dechin Chen, Zhijun Pan, Yu Li, Zijing Liu, Yi Isaac Yang.   
[arXiv:2504.18367 (2025)](https://arxiv.org/abs/2504.18367) | [data](https://huggingface.co/SZBL-IDEA)  

* **QMe14S: A Comprehensive and Efficient Spectral Data Set for Small Organic Molecules** [2025]  
Cristian Gabellini, Nikhil Shenoy, Stephan Thaler, Semih Canturk, Daniel McNeela, Dominique Beaini, Michael Bronstein, Prudencio Tossou.   
[J. Phys. Chem. Lett. (2025)](https://doi.org/10.1021/acs.jpclett.5c00839) | [data](https://figshare.com/s/889262a4e999b5c9a5b3)  

* **UniSim: A Unified Simulator for Time-Coarsened Dynamics of Biomolecules** [2025]  
Ziyang Yu, Wenbing Huang, Yang Liu.  
[ICML 2025 (2025)](https://doi.org/10.48550/arXiv.2506.03157) | [code&data](https://zenodo.org/records/15394209)  

* **The QCML dataset, Quantum chemistry reference data from 33.5M DFT and 14.7B semi-empirical calculations** [2025]  
Ganscha, S., Unke, O.T., Ahlin, D. et al.  
[Sci Data 12, 406 (2025)](https://doi.org/10.1038/s41597-025-04720-7) | [data](https://zenodo.org/records/14859804)  

* **OpenQDC: Open Quantum Data Commons** [2024]  
Cristian Gabellini, Nikhil Shenoy, Stephan Thaler, Semih Canturk, Daniel McNeela, Dominique Beaini, Michael Bronstein, Prudencio Tossou.   
[arXiv:2411.19629 (2024)](https://arxiv.org/abs/2411.19629) | [data](https://github.com/valence-labs/openQDC)  

* **Molecular Quantum Chemical Data Sets and Databases for Machine Learning Potentials** [2024]  
Antonio Mirarchi, Toni Giorgino, G. D. Fabritiis.   
[ChemRxiv. (2024)](https://doi.org/10.26434/chemrxiv-2024-w3ld0-v2) | [code](https://github.com/Arif-PhyChem/datasets_and_databases_4_MLPs)  

* **mdCATH: A Large-Scale MD Dataset for Data-Driven Computational Biophysics** [2024]  
Antonio Mirarchi, Toni Giorgino, G. D. Fabritiis.   
[	arXiv:2407.14794 (2024)](https://arxiv.org/abs/2407.14794) | [code](https://github.com/compsciencelab/mdCATH)  

* **nablaDFT: Large-Scale Conformational Energy and Hamiltonian Prediction benchmark and dataset** [2022]  
Khrabrov, Kuzma and Shenbin, Ilya and Ryabov, Alexander and Tsypin, Artem and Telepov, Alexander and Alekseev, Anton and Grishin, Alexander and Strashnov, Pavel and Zhilyaev, Petr and Nikolenko, Sergey and Kadurin, Artur.   
[Phys. Chem. Chem. Phys. (2022)](https://doi.org/10.1039/D2CP03966D) | [code](https://github.com/AIRI-Institute/nablaDFT)  




### Package




**MMolearn**  
a Python package streamlining the design of generative models of biomolecular dynamics  

https://github.com/LumosBio/MolData   







## Molecular dynamics



### Molecular Force Fields




* **ABEEM Polarizable Force Field for PC Lipids: Parameterization and Molecular Dynamics Simulations** [2025]  
Xiaoyu Wang, Linlin Liu, Peiran Meng, Jian Zhao, Lei Wang, Cui Liu, Lidong Gong, and Zhongzhi Yang.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.4c01704)  

* **The PHAST 2.0 Force Field for General Small Molecule and Materials Simulations** [2025]  
Adam Hogan, Logan Ritter, and Brian Space.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.5c00134)  

* **Martini3-IDP: improved Martini 3 force field for disordered proteins** [2025]  
Wang, L., Brasnett, C., Borges-Araújo, L. et al.  
[Nat Commun 16, 2874 (2025)](https://doi.org/10.1038/s41467-025-58199-2) | [code](https://github.com/Martini-Force-Field-Initiative/Martini3-IDP-parameters)  






### Neural Molecular Force Fields








* **aims-PAX: Parallel Active eXploration for the automated construction of Machine Learning Force Fields** [2025]  
Tobias Henkes, Shubham Sharma, Alexandre Tkatchenko, Mariana Rossi, Igor Poltavskyi.  
[arXiv:2508.12888 (2025)](https://doi.org/10.48550/arXiv.2508.12888) | [code](https://github.com/tohenkes/aims-PAX)  

* **A density based machine learning force field for molecule-molecule nonbonded interactions** [2025]  
Wang L-W.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-83grg)  

* **Force field optimization by end-to-end differentiable atomistic simulation** [2025]  
Abhijeet Sadashiv Gangan, Ekin Dogus Cubuk, Samuel S. Schoenholz, Mathieu Bauchy, and N. M. Anoop Krishnan.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.4c01784) | [code](https://github.com/M3RG-IITD/Force-field-optimization)  

* **Operator Forces For Coarse-Grained Molecular Dynamics** [2025]  
Klein, Leon, Atharva Kelkar, Aleksander Durumeric, Yaoyi Chen, and Frank Noé.  
[arXiv:2506.19628 (2025)](https://doi.org/10.48550/arXiv.2506.19628) | [code](https://github.com/noegroup/OperatorForces4CG)  

* **Evolutionary machine learning of physics-based force fields in high-dimensional parameter-space** [2025]  
van der Spoel D, Marrades J, Kriz K, Hosseini AN, Nordman A, Ateide Martins JP, et al.  
[	Digital Discovery (2025)](https://doi.org/10.1039/D5DD00178A) | [code](https://github.com/AlexandriaChemistry/ACT)  

* **To Use or Not to Use a Universal Force Field** [2025]  
Denan Li, Jiyuan Yang, Xiangkai Chen, Lintao Yu, Shi Liu.  
[arXiv:2503.08207 (2025)](https://arxiv.org/abs/2503.08207)  

* **Understanding and Mitigating Distribution Shifts For Machine Learning Force Fields** [2025]  
Kreiman, Tobias, and Aditi S. Krishnapriyan.  
[arXiv:2503.08674 (2025)](https://arxiv.org/abs/2503.08674) | [code](https://github.com/ASK-Berkeley/MLFF-distribution-shifts)  

* **Accelerating CO2 direct air capture screening for metal-organic frameworks with a transferable machine learning force field** [2025]  
Yunsung Lim and Hyunsoo Park and Aron Walsh and Jihan Kim.  
[Matter (2025)](https://doi.org/10.1016/j.matt.2025.102203) | [code](https://github.com/hspark1212/DAC-SIM)  

* **On the design space between molecular mechanics and machine learning force fields** [2025]  
Wang, Yuanqing, Kenichiro Takaba, Michael S. Chen, Marcus Wieder, Yuzhi Xu, Tong Zhu, John ZH Zhang et al.  
[Appl. Phys. Rev. (2025)](https://doi.org/10.1063/5.0237876)  

* **ABFML: A problem-oriented package for rapidly creating, screening, and optimizing new machine learning force fields** [2025]  
Xingze Geng, Jianing Gu, Gaowu Qin, Lin-Wang Wang, Xiangying Meng.  
[J. Chem. Phys. (2025)](https://doi.org/10.1063/5.0247559) | [code](https://github.com/gengxingze/ABFML)  

* **Reversible molecular simulation for training classical and machine-learning force fields** [2025]  
Greener, Joe G.   
[Proc. Natl. Acad. Sci. (2025)](https://doi.org/10.1073/pnas.2426058122) | [code](https://github.com/greener-group/rev-sim)  

* **Artificial Intelligence for Direct Prediction of Molecular Dynamics Across Chemical Space** [2025]  
Ge F, Dral PO.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-kc7sn) | [code](https://github.com/dralgroup/mlatom)  

* **NepoIP/MM: Toward Accurate Biomolecular Simulation with a Machine Learning/Molecular Mechanics Model Incorporating Polarization Effects** [2025]  
Ge Song and Weitao Yang.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.5c00372) | [code](https://github.com/Yang-Laboratory/NepoIP)  

* **Efficient Long-Range Machine Learning Force Fields for Liquid and Materials Properties** [2025]  
Weber, John L., Rishabh D. Guha, Garvit Agarwal, Yujing Wei, Aidan A. Fike, Xiaowei Xie, James Stevenson et al.  
[arXiv:2505.06462 (2025)](https://arxiv.org/abs/2505.06462) | [code](https://github.com/leifjacobson/MLFF_test_data)  

* **MACE-OFF: Short-Range Transferable Machine Learning Force Fields for Organic Molecules** [2025]  
Dávid Péter Kovács, J. Harry Moore, Nicholas J. Browning, Ilyes Batatia, Joshua T. Horton, Yixuan Pu, Venkat Kapil, William C. Witt, Ioan-Bogdan Magdău, Daniel J. Cole, and Gábor Csányi.  
[J. Am. Chem. Soc. (2025)](https://doi.org/10.1021/jacs.4c07099) | [code](https://github.com/ACEsuit/mace-off)  

* **Grappa--A Machine Learned Molecular Mechanics Force Field** [2025]  
Seute, Leif, Eric Hartmann, Jan Stühmer, and Frauke Gräter.  
[Chemical Science 16.6 (2025)](https://doi.org/10.1039/D4SC05465B) | [arXiv:2404.00050 (2024)](https://arxiv.org/abs/2404.00050) | [code](https://github.com/graeter-group/grappa)  

* **ILVES: Accurate and efficient bond length and angle constraints in molecular dynamics** [2025]  
López-Villellas, L., Mikkelsen, C.C.K., Galano-Frutos, J.J., Marco-Sola, S., Alastruey-Benedé, J., Ibáñez, P., Moretó, M., De Rosa, M.C. and García-Risueño, P.  
[arXiv:2503.13075 (2025)](https://arxiv.org/abs/2503.13075)  

* **Towards Fast, Specialized Machine Learning Force Fields: Distilling Foundation Models via Energy Hessians** [2025]  
Ishan Amin, Sanjeev Raja, Aditi Krishnapriyan.  
[arXiv:2501.09009 (2025)](https://arxiv.org/abs/2501.09009) | [code](https://github.com/ASK-Berkeley/MLFF-distill)  

* **BoostMD: Accelerated Molecular Sampling Leveraging ML Force Field Features** [2025]  
Zhaoxin Xie, Yanheng Li, Yijie Xia, Jun Zhang, Sihao Yuan, Cheng Fan, Yi Isaac Yang, and Yi Qin Gao.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.4c01449)  

* **Efficient and Precise Force Field Optimization for Biomolecules Using DPA-2** [2024]  
Junhan Chang, Duo Zhang, Yuqing Deng, Hongrui Lin, Zhirong Liu, Linfeng Zhang, Hang Zheng, Xinyan Wang.   
[arXiv:2406.09817 (2024)](https://arxiv.org/abs/2406.09817) | [code](https://github.com/deepmodeling/DeePDih)  

* **FeNNol: An efficient and flexible library for building force-field-enhanced neural network potentials** [2024]  
Thomas Plé, Olivier Adjoua, Louis Lagardère, Jean-Philip Piquemal.   
[J. Chem. Phys. (2024)](https://doi.org/10.1063/5.0217688) | [code](https://github.com/thomasple/FeNNol)  

* **Reversible molecular simulation for training classical and machine learning force fields** [2024]  
Joe G Greener.   
[arXiv:2412.04374 (2024)](https://arxiv.org/abs/2412.04374) | [code](https://github.com/JuliaMolSim/Molly.jl)  

* **HessFit: A Toolkit to Derive Automated Force Fields from Quantum Mechanical Information** [2024]  
Falbo, E. and Lavecchia, A.   
[J. Chem. Inf. Model. (2024)](https://doi.org/10.1021/acs.jcim.4c00540) | [code](https://github.com/emanuelefalbo/HessFit)  

* **A Euclidean transformer for fast and stable machine learned force fields** [2024]  
Frank, J.T., Unke, O.T., Müller, KR. et al.   
[Nat Commun 15, 6539 (2024)](https://doi.org/10.1038/s41467-024-50620-6) | [code](https://github.com/microsoft/AI2BMD/tree/ViSNet/chignolin_data)

* **Differentiable simulation to develop molecular dynamics force fields for disordered proteins** [2024]  
Greener, Joe G.   
[Chemical Science 15.13 (2024)](https://doi.org/10.1039/D3SC05230C) | [code](https://github.com/greener-group/GB99dms)  

* **Grappa--A Machine Learned Molecular Mechanics Force Field** [2024]  
Seute, Leif, Eric Hartmann, Jan Stühmer, and Frauke Gräter.  
[arXiv:2404.00050 (2024)](https://arxiv.org/abs/2404.00050) | [code](https://github.com/graeter-group/grappa)  

* **An implementation of the Martini coarse-grained force field in OpenMM** [2023]  
MacCallum, J. L., Hu, S., Lenz, S., Souza, P. C., Corradi, V., & Tieleman, D. P.   
[Biophysical Journal 122.14 (2023)](https://doi.org/10.1016/j.bpj.2023.04.007)  

* **Improving machine learning force fields for molecular dynamics simulations with fine-grained force metrics** [2023]  
Zun Wang, Hongfei Wu, Lixin Sun, Xinheng He, Zhirong Liu, Bin Shao, Tong Wang, Tie-Yan Liu.   
[J. Chem. Phys. (2023)](https://doi.org/10.1063/5.0147023) | [data](https://figshare.com/articles/dataset/Revised_MD17_dataset_rMD17_/12672038)    

* **End-to-end differentiable construction of molecular mechanics force fields** [2022]  
Wang, Yuanqing, Josh Fass, Benjamin Kaminow, John E. Herr, Dominic Rufa, Ivy Zhang, Iván Pulido et al.   
[Chemical Science (2022)](https://doi.org/10.1039/D2SC02739A) | [code](https://github.com/openforcefield/openff-nagl)  

* **Improving molecular force fields across configurational space by combining supervised and unsupervised machine learning** [2021]  
Gregory Fonseca, Igor Poltavsky, Valentin Vassilev-Galindo, Alexandre Tkatchenko.   
[J. Chem. Phys. (2021)](https://doi.org/10.1063/5.0035530)  




### MD Engines-Frameworks

* [Amber](http://ambermd.org/) - A suite of biomolecular simulation programs.  
* [Gromacs](http://www.gromacs.org/) - A molecular dynamics package mainly designed for simulations of proteins, lipids and nucleic acids.
* [OpenMM](http://openmm.org/) - A toolkit for molecular simulation using high performance GPU code.  
* [CHARMM](https://www.charmm.org/) - A molecular simulation program with broad application to many-particle systems.  
* [HTMD](https://github.com/Acellera/htmd) - Programming Environment for Molecular Discovery.  
* [ACEMD](https://www.acellera.com/acemd) - The next generation molecular dynamic simulation software.  
* [NAMD](https://www.ks.uiuc.edu/Research/namd/) - A parallel molecular dynamics code for large biomolecular systems.
* [StreaMD](https://github.com/ci-lab-cz/streamd) - A tool to perform high-throughput automated molecular dynamics simulations.
* [BEMM-GEN](https://github.com/y4suda/BEMM-GEN) - A Toolkit for Generating a Biomolecular Environment-Mimicking Model for Molecular Dynamics Simulation.
* [BioSimSpace](https://github.com/OpenBioSim/biosimspace) - An interoperable Python framework for biomolecular simulation.
* [STORMM](https://github.com/Psivant/stormm) - Structure and TOpology Replica Molecular Mechanics.




* **OpenCafeMol: A GPU-accelerated coarse-grained biomolecular dynamics simulator with OpenMM library** [2025]  
Barnett, Simon, and John D. Chodera.  
[bioRxiv. (2025)](https://doi.org/10.1101/2025.02.20.639390)  





### AI4MD Engines-Frameworks

* [OpenMM 8](https://github.com/openmm/openmm) - Molecular Dynamics Simulation with Machine Learning Potentials.  
* [DeePMD-kit](https://github.com/deepmodeling/deepmd-kit) - A deep learning package for many-body potential energy representation and molecular dynamics.  
* [TorchMD](https://github.com/torchmd/torchmd) - End-To-End Molecular Dynamics (MD) Engine using PyTorch.  
* [TorchMD-NET](https://github.com/torchmd/torchmd-net) - TorchMD-NET provides state-of-the-art neural networks potentials (NNPs) and a mechanism to train them.  
* [OpenMM-Torch](https://github.com/openmm/openmm-torch) - OpenMM plugin to define forces with neural networks.  
* [AI<sup>2</sup>BMD](https://github.com/microsoft/AI2BMD) - AI-powered ab initio biomolecular dynamics simulation.  
* [NeuralMD](https://github.com/chao1224/NeuralMD) - A Multi-Grained Symmetric Differential Equation Model for Learning Protein-Ligand Binding Dynamics.  
* [TorchSim](https://github.com/Radical-AI/torch-sim) - A next-generation open-source atomistic simulation engine for the MLIP era.  




### MD Trajectory Processing-Analysis


* [MDAnalysis](https://www.mdanalysis.org/) - An object-oriented Python library to analyze trajectories from molecular dynamics (MD) simulations in many popular formats.  
* [MDTraj](http://mdtraj.org/) - A python library that allows users to manipulate molecular dynamics (MD) trajectories.  
* [PyTraj](https://amber-md.github.io/pytraj/) - A Python front-end package of the popular cpptraj program.  
* [CppTraj](https://github.com/Amber-MD/cpptraj) - Biomolecular simulation trajectory/data analysis.  
* [WEDAP](https://github.com/chonglab-pitt/wedap) - A Python Package for Streamlined Plotting of Molecular Simulation Data. 
* [Melodia](https://github.com/rwmontalvao/Melodia_py) - A Python library for protein structure analysis. 
* [MDANCE](https://github.com/mqcomplab/MDANCE) - A flexible n-ary clustering package that provides a set of tools for clustering Molecular Dynamics trajectories.  
* [PENSA](https://github.com/drorlab/pensa) - A collection of python methods for exploratory analysis and comparison of biomolecular conformational ensembles.  






* **apoCHARMM: High-performance molecular dynamics simulations on GPUs for advanced simulation methods** [2025]  
Samarjeet Prasad, Felix Aviat, James E. Gonzales, Bernard R. Brooks.  
[J. Chem. Phys. (2025)](https://doi.org/10.1063/5.0264937) | [code](https://github.com/samarjeet/apocharmm)  

* **Self-Supervised Evolution Operator Learning for High-Dimensional Dynamical Systems** [2025]  
Giacomo Turri, Luigi Bonati, Kai Zhu, Massimiliano Pontil, Pietro Novelli.  
[arXiv:2505.18671 (2025)](https://doi.org/10.48550/arXiv.2505.18671) | [code](https://github.com/pietronvll/encoderops)  

* **NeuralTSNE: A Python Package for the Dimensionality Reduction of Molecular Dynamics Data Using Neural Networks** [2025]  
Patryk Tajs, Mateusz Skarupski, Jakub Rydzewski.  
[arXiv:2505.16476(2025)](https://doi.org/10.48550/arXiv.2505.16476) | [code](https://github.com/NeuralTSNE/NeuralTSNE)  

* **GEODES: Geometric Descriptors for the Assessment of Global and Local Flexibility of Proteins During Molecular Dynamics Simulation** [2025]  
Pats, Karina and Glukhov, Igor and Petrosian, Stepan and Mamaeva, Maria and Sergushichev, Alexey and Devignes, Marie-Dominique and Molnár, Ferdinand.  
[IEEE Access (2025)](https://doi.org/10.1109/ACCESS.2025.3558781) | [code](https://github.com/rinnifox/GEODES)  

* **gmx_RRCS: a precision tool for detecting subtle conformational dynamics in molecular simulations** [2025]  
Wei Han, Zhenghan Chen, Ming-Wei Wang, Qingtong Zhou.  
[Journal of Molecular Biology (2025)](https://doi.org/10.1016/j.jmb.2025.169129) | [code](https://github.com/RuijinHospitalRCMSB/gmx_RRCS)  

* **Systematic analysis of biomolecular conformational ensembles with PENSA** [2025]  
Vögele, Martin, Neil J. Thomson, Sang T. Truong, Jasper McAvity, Ulrich Zachariae, and Ron O. Dror.  
[The Journal of Chemical Physics 162.1 (2025)](https://doi.org/10.1063/5.0235544) | [code](https://github.com/drorlab/pensa)  

* **MDRefine: a Python package for refining Molecular Dynamics trajectories with experimental data** [2024]  
Ivan Gilardoni, Valerio Piomponi, Thorben Fröhlking, Giovanni Bussi.   
[	arXiv:2411.07798 (2024)](https://arxiv.org/abs/2411.07798) | [code](https://github.com/bussilab/MDRefine)  

* **NRIMD, a Web Server for Analyzing Protein Allosteric Interactions Based on Molecular Dynamics Simulation** [2024]  
He, Yi, Shuang Wang, Shuai Zeng, Jingxuan Zhu, Dong Xu, Weiwei Han, and Juexin Wang.   
[J. Chem. Inf. Model. (2024)](https://doi.org/10.1021/acs.jcim.4c00783) | [web](https://nrimd.luddy.indianapolis.iu.edu/)  





### Reference

https://github.com/ipudu/awesome-molecular-dynamics  



### Visualization

* [VMD](http://www.ks.uiuc.edu/Research/vmd/) - A molecular visualization program for displaying, animating, and analyzing large biomolecular systems using 3-D graphics and built-in scripting.  
* [NGLview](https://github.com/arose/nglview) - IPython widget to interactively view molecular structures and trajectories.  
* [PyMOL](https://pymol.org/2/) - A user-sponsored molecular visualization system on an open-source foundation, maintained and distributed by Schrödinger.  
* [Avogadro](https://avogadro.cc/) - An advanced molecule editor and visualizer designed for cross-platform use in computational chemistry, molecular modeling, bioinformatics, materials science, and related areas.  






## CGMD
**Coarse-Grained Molecular Dynamics Simulations**






* **Understanding Viscoelasticity of an Entangled Silicone Copolymer via Coarse-Grained Molecular Dynamics Simulations** [2025]  
Weikang Xian, Amitesh Maiti, Andrew P. Saab, and Ying Li.  
[Macromolecules (2025)](https://doi.org/10.1021/acs.macromol.5c01192)  

* **Navigating protein landscapes with a machine-learned transferable coarse-grained model** [2025]  
Charron, N.E., Bonneau, K., Pasos-Trejo, A.S. et al.  
[Nat. Chem. (2025)](https://doi.org/10.1038/s41557-025-01874-0) | [data](https://github.com/ClementiGroup/mlcg)  

* **Graph-Coarsening for Machine Learning Coarse-grained Molecular Dynamics** [2025]  
Soumya Mondal, Subhanu Halder, Debarchan Basu, Sandeep Kumar, Tarak Karmakar.  
[	arXiv:2507.16531 (2025)](https://doi.org/10.48550/arXiv.2507.16531)  

* **Fast parameterization of Martini3 models for fragments and small molecules** [2025]  
Magdalena Szczuka, Gilberto P. Pereira, Luis J. Walter, Marc Gueroult, Pierre Poulain, Tristan Bereau, Paulo C. T. Souza, Matthieu Chavent.  
[bioRxiv (2025)](https://doi.org/10.1101/2025.07.13.664596) | [data](https://github.com/Martini-Force-Field-Initiative/Automartini_M3)  






## AI4MD


* **TorchSim: An efficient atomistic simulation engine in PyTorch** [2025]  
Orion Cohen, Janosh Riebesell, Rhys Goodall, Adeesh Kolluru, Stefano Falletta, Joseph Krause, Jorge Colindres, Gerbrand Ceder, Abhijeet S. Gangan.  
[arXiv:2508.06628 (2025)](https://doi.org/10.48550/arXiv.2508.06628) | [code](https://github.com/Radical-AI/torch-sim)  

* **chemtrain-deploy: A parallel and scalable framework for machine learning potentials in million-atom MD simulations** [2025]  
Paul Fuchs, Weilong Chen, Stephan Thaler, Julija Zavadlav.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.5c00996) | [arXiv:2506.04055 (2025)](https://doi.org/10.48550/arXiv.2506.04055) | [data](https://github.com/tummfm/chemtrain)  

* **Action-Minimization Meets Generative Modeling: Efficient Transition Path Sampling with the Onsager-Machlup Functional** [2025]  
Sanjeev Raja, Martin Šípka, Michael Psenka, Tobias Kreiman, Michal Pavelka, Aditi S. Krishnapriyan.  
[ICML 2025 (2025)](https://arxiv.org/abs/2504.18506) | [data](https://github.com/ASK-Berkeley/OM-TPS)  

* **Operator Forces For Coarse-Grained Molecular Dynamics** [2025]  
Leon Klein, Atharva Kelkar, Aleksander Durumeric, Yaoyi Chen, Frank Noé.  
[arXiv:2506.19628 (2025)](https://doi.org/10.48550/arXiv.2506.19628) | [data](https://github.com/noegroup/OperatorForces4CG)  

* **Memory kernel minimization-based neural networks for discovering slow collective variables of biomolecular dynamics** [2025]  
Liu, B., Cao, S., Boysen, J.G. et al.  
[Nat Comput Sci (2025)](https://doi.org/10.1038/s43588-025-00815-8) | [code](https://github.com/markovmodel/mdshare)  

* **UniSim: A Unified Simulator for Time-Coarsened Dynamics of Biomolecules** [2025]  
Ziyang Yu, Wenbing Huang, Yang Liu.  
[ICML 2025 (2025)](https://doi.org/10.48550/arXiv.2506.03157) | [code](https://github.com/yaledeus/UniSim) | [data](https://zenodo.org/records/15394209)  

* **FlashMD: long-stride, universal prediction of molecular dynamics** [2025]  
Bigi, Filippo, Sanggyu Chong, Agustinus Kristiadi, and Michele Ceriotti.  
[arXiv:2505.19350 (2025)](https://arxiv.org/abs/2505.19350) | [code](https://github.com/lab-cosmo/flashmd)  

* **Investigating the Nature of PRM:SH3 Interactions Using Artificial Intelligence and Molecular Dynamics** [2025]  
Se-Jun Kim, Da-Eun Hwang, Hyungjun Kim, and Jeong-Mo Choi.  
[J. Chem. Inf. Model. (2025)](https://doi.org/10.1021/acs.jcim.5c00342) | [code](https://github.com/sejunkim6370/AlphaFold2-and-Molecular-Dynamics)  

* **Towards Unraveling Biomolecular Conformational Landscapes with a Generative Foundation Model** [2025]  
The OpenComplex Team, Qiwei Ye.  
[bioRxiv. (2025)](https://doi.org/10.1101/2025.05.01.651643)  

* **Deep Signature: Characterization of Large-Scale Molecular Dynamics** [2025]  
Tiexin Qin, Mengxu ZHU, Chunyang Li, Terry Lyons, Hong Yan, Haoliang Li.  
[ICLR (2025)](https://openreview.net/forum?id=xayT1nn8Mg) | [code](https://github.com/WonderSeven/Deep-Signature)  

* **A Foundation Model for Accurate Atomistic Simulations in Drug Design** [2025]  
Plé T, Adjoua O, Benali A, Posenitskiy E, Villot C, Lagardère L, et al.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-f1hgn) | [code](https://github.com/thomasple/FeNNol-PMC)  

* **Orb-v3: atomistic simulation at scale** [2025]  
Benjamin Rhodes, Sander Vandenhaute, Vaidotas Šimkus, James Gin, Jonathan Godwin, Tim Duignan, Mark Neumann.  
[arXiv:2504.06231 (2025)](https://arxiv.org/abs/2504.06231) | [code](https://github.com/orbital-materials/orb-models)  

* **A predictive machine learning force-field framework for liquid electrolyte development** [2025]  
Gong, S., Zhang, Y., Mu, Z. et al.  
[Nat Mach Intell (2025)](https://doi.org/10.1038/s42256-025-01009-7) | [code](https://github.com/bytedance/bamboo)  

* **Fast, Modular, and Differentiable Framework for Machine Learning-Enhanced Molecular Simulations** [2025]  
Henrik Christiansen, Takashi Maruyama, Federico Errica, Viktor Zaverkin, Makoto Takamoto, Francesco Alesiani.  
[arXiv:2503.20541 (2025)](https://arxiv.org/abs/2503.20541) | [code](https://github.com/nec-research/DIMOS)  

* **Kinetically Consistent Coarse Graining using Kernel-based Extended Dynamic Mode Decomposition** [2025]  
Vahid Nateghi, Feliks Nüske.   
[arXiv:2409.16396 (2025)](https://arxiv.org/abs/2409.16396) | [code](https://zenodo.org/records/13808938)  

* **Advancing nonadiabatic molecular dynamics simulations in solids with E(3) equivariant deep neural hamiltonians** [2025]  
Zhang, C., Zhong, Y., Tao, ZG. et al.   
[Nat Commun 16, 2033 (2025)](https://doi.org/10.1038/s41467-025-57328-1) | [code](https://github.com/QuantumLab-ZY/HamGNN)  

* **BoostMD: Accelerated Molecular Sampling Leveraging ML Force Field Features** [2024]  
Schaaf, Lars Leon, Ilyes Batatia, Jules Tilly, and Thomas D. Barrett.   
[NeurIPS 2024 Workshop on Data-driven and Differentiable Simulations, Surrogates, and Solvers. (2024)](https://openreview.net/forum?id=H0USH61HnF)  

* **A Multi-Grained Symmetric Differential Equation Model for Learning Protein-Ligand Binding Dynamics** [2024]  
Shengchao Liu, Weitao Du, Hannan Xu, Yanjing Li, Zhuoxinran Li, Vignesh Bhethanabotla, Divin Yan, Christian Borgs, Anima Anandkumar, Hongyu Guo, Jennifer Chayes.   
[arXiv:2401.15122 (2024)](https://arxiv.org/abs/2401.15122) | [code](https://github.com/chao1224/NeuralMD)  

* **Accelerating Molecular Dynamics Simulations with Quantum Accuracy by Hierarchical Classification** [2024]  
Brunken, Christoph, Sebastien Boyer, Mustafa Omar, Bakary N'tji Diallo, Karim Beguir, Nicolas Lopez Carranza, and Oliver Bent.   
[ChemRxiv. (2024)](https://doi.org/10.26434/chemrxiv-2024-20q1v) | [code](https://github.com/multiscale-modelling/humf)  

* **Machine learning of force fields towards molecular dynamics simulations of proteins at DFT accuracy** [2024]  
Brunken, Christoph, Sebastien Boyer, Mustafa Omar, Bakary N'tji Diallo, Karim Beguir, Nicolas Lopez Carranza, and Oliver Bent.   
[ICLR 2024 Workshop on Generative and Experimental Perspectives for Biomolecular Design (2024)](https://openreview.net/forum?id=hrvvIOx7EM) | [code](https://github.com/ACEsuit/mace-jax)  

* **Unsupervised and supervised AI on molecular dynamics simulations reveals complex characteristics of HLA-A2-peptide immunogenicity** [2024]  
Jeffrey K Weber, Joseph A Morrone, Seung-gu Kang, Leili Zhang, Lijun Lang, Diego Chowell, Chirag Krishna, Tien Huynh, Prerana Parthasarathy, Binquan Luan, Tyler J Alban, Wendy D Cornell, Timothy A Chan.   
[Briefings in Bioinformatics (2024)](https://doi.org/10.1093/bib/bbad504) | [code](https://github.com/BiomedSciAI)

* **Biomolecular dynamics with machine-learned quantum-mechanical force fields trained on diverse chemical fragments** [2024]  
Unke, Oliver T., Martin Stöhr, Stefan Ganscha, Thomas Unterthiner, Hartmut Maennel, Sergii Kashubin, Daniel Ahlin et al.   
[Science Advances 10.14 (2024)](https://www.science.org/doi/10.1126/sciadv.adn4397) | [data](https://zenodo.org/records/10720941)

* **DeePMD-kit v2: A software package for deep potential models** [2023]  
Zeng, Jinzhe, Duo Zhang, Denghui Lu, Pinghui Mo, Zeyu Li, Yixiao Chen, Marián Rynik et al.   
[The Journal of Chemical Physics 159.5 (2023)](https://doi.org/10.1063/5.0155600) | [code](https://github.com/deepmodeling/deepmd-kit)  

* **Enhancing Biomolecular Sampling with Reinforcement Learning: A Tree Search Molecular Dynamics Simulation Method** [2019]  
Kento Shin, Duy Phuoc Tran, Kazuhiro Takemura, Akio Kitao, Kei Terayama, and Koji Tsuda.   
[ACS Omega (2019)](https://doi.org/10.1021/acsomega.9b01480) | [code](https://github.com/tsudalab/TSMD)  

* **DeePMD-kit: A deep learning package for many-body potential energy representation and molecular dynamics** [2018]  
Wang, Han, Linfeng Zhang, Jiequn Han, and E. Weinan.   
[Computer Physics Communications 228 (2018)](https://doi.org/10.1016/j.cpc.2018.03.016) | [code](https://github.com/deepmodeling/deepmd-kit)  











### Neural Network Potentials







* **A critical review of machine learning interatomic potentials and Hamiltonian** [2025]  
Li, Y.; Zhang, X.; Liu, M.; Shen, L.  
  [J. Mater. Inf. (2025)](http://dx.doi.org/10.20517/jmi.2025.17)  

* **PIL-Net: a physics-informed graph convolutional network for predicting atomic multipoles** [2025]  
Dario Coscia, Pim de Haan, Max Welling.  
[arXiv:2508.14022 (2025)](https://arxiv.org/abs/2508.14022) | [code](https://github.com/dario-coscia/blip)  

* **FastTrack: a fast method to evaluate mass transport in solid leveraging universal machine learning interatomic potential** [2025]  
Hanwen Kang, Tenglong Lu, Zhanbin Qi, Jiandong Guo, Sheng Meng, Miao Liu.  
[arXiv:2508.10505 (2025)](https://doi.org/10.48550/arXiv.2508.10505)  

* **PIL-Net: a physics-informed graph convolutional network for predicting atomic multipoles** [2025]  
Whitter, Caitlin, Alex Pothen, and Aurora Clark.  
[Digital Discovery (2025)](https://doi.org/10.1039/D5DD00228A) | [code](https://github.com/rinikerlab/EquivariantMultipoleGNN)  

* **Performance of universal machine-learned potentials with explicit long-range interactions in biomolecular simulations** [2025]  
Viktor Zaverkin, Matheus Ferraz, Francesco Alesiani, Mathias Niepert.  
[arXiv:2508.10841 (2025)](https://doi.org/10.48550/arXiv.2508.10841) | [code](https://github.com/nec-research/ictp)  

* **TorchANI 2.0: An extensible, high performance library for the design, training, and use of NN-IPs** [2025]  
Pickering I, Xue J, Huddleston K, Terrel N, Roitberg A.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-w0w7d) | [code](https://github.com/aiqm/torchani)  

* **chemtrain-deploy: A parallel and scalable framework for machine learning potentials in million-atom MD simulations** [2025]  
Paul Fuchs, Weilong Chen, Stephan Thaler, Julija Zavadlav.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.5c00996) | [arXiv:2506.04055 (2025)](https://doi.org/10.48550/arXiv.2506.04055) | [data](https://github.com/tummfm/chemtrain)  

* **Augmenting Molecular Graphs with Geometries via Machine Learning Interatomic Potentials** [2025]  
Fu, Cong, Yuchao Lin, Zachary Krueger, Haiyang Yu, Maho Nakata, Jianwen Xie, Emine Kucukbenli, Xiaofeng Qian, and Shuiwang Ji.  
[arXiv:2507.00407 (2025)](https://doi.org/10.48550/arXiv.2507.00407) | [data](https://huggingface.co/divelab)  

* **Tensor Decomposition Networks for Fast Machine Learning Interatomic Potential Computations** [2025]  
Lin, Yuchao, Cong Fu, Zachary Krueger, Haiyang Yu, Maho Nakata, Jianwen Xie, Emine Kucukbenli, Xiaofeng Qian, and Shuiwang Ji.  
[arXiv:2507.01131 (2025)](https://doi.org/10.48550/arXiv.2507.01131) | [data](https://huggingface.co/divelab)  

* **A Benchmark for Quantum Chemistry Relaxations via Machine Learning Interatomic Potentials** [2025]  
Fu, Cong, Yuchao Lin, Zachary Krueger, Wendi Yu, Xiaoning Qian, Byung-Jun Yoon, Raymundo Arróyave et al.  
[arXiv:2506.23008 (2025)](https://doi.org/10.48550/arXiv.2506.23008) | [data](https://huggingface.co/divelab)  

* **Enhancing Machine Learning Potentials through Transfer Learning across Chemical Elements** [2025]  
Sebastien Röcken and Julija Zavadlav.  
[J. Chem. Inf. Model. (2025)](https://doi.org/10.1021/acs.jcim.5c00293) | [code](https://github.com/tummfm/TL_MLP)  

* **Distillation of atomistic foundation models across architectures and chemical domains** [2025]  
Gardner, John LA, Daniel F. Toit, Chiheb Ben Mahmoud, Zoé Faure Beaulieu, Veronika Juraskova, Laura-Bianca Paşca, Louise AM Rosset et al.  
[arXiv:2506.10956 (2025)](https://arxiv.org/abs/2506.10956) | [code](https://github.com/jla-gardner/augment-atoms)  

* **Spin-informed universal graph neural networks for simulating magnetic ordering** [2025]  
W. Xu,R.Y. Sanspeur,A. Kolluru,B. Deng,P. Harrington,S. Farrell,K. Reuter,& J.R. Kitchin.  
[Proc. Natl. Acad. Sci. (2025)](https://doi.org/10.1073/pnas.2422973122) | [code](https://github.com/Wenbintum/ocp_mag)  

* **UMA: A Family of Universal Models for Atoms** [2025]  
Brandon M. Wood, Misko Dzamba, Xiang Fu, Meng Gao, Muhammed Shuaibi, Luis Barroso-Luque, Kareem Abdelmaqsoud, Vahe Gharakhanyan, John R. Kitchin, Daniel S. Levine, Kyle Michel, Anuroop Sriram, Taco Cohen, Abhishek Das, Ammar Rizvi, Sushree Jagriti Sahoo, Zachary W. Ulissi, C. Lawrence Zitnick.  
[arXiv:2506.23971(2025)](https://doi.org/10.48550/arXiv.2506.23971) | [code](https://huggingface.co/facebook/UMA)  

* **Application-specific Machine-Learned Interatomic Potentials: Exploring the Trade-off Between Precision and Computational Cost** [2025]  
Baghishov, Ilgar, Jan Janssen, Graeme Henkelman, and Danny Perez.  
[arXiv:2506.05646 (2025)](https://arxiv.org/abs/2506.05646)  

* **DeePEST-OS: A Generic Machine Learning Potential for Accelerating Transition State Search in Organic Synthesis** [2025]  
Ren K, Tang K, Zhao Y, Zhang L, Du J, Meng Q, et al.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-mzz6w)  

* **AIMNet2-rxn: A Machine Learned Potential for Generalized Reaction Modeling on a Millions-of-Pathways Scale** [2025]  
Anstine DM, Zhao Q, Zubatiuk R, Zhang S, Singla V, Nikitin F, et al.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-hpdmg) | [code](https://github.com/isayevlab/aimnetcentral)  

* **Global universal scaling and ultrasmall parameterization in machine-learning interatomic potentials with superlinearity** [2025]  
Y. Hu,Y. Sheng,J. Huang,X. Xu,Y. Yang,M. Zhang,Y. Wu,C. Ye,J. Yang,& W. Zhang,  
[Proc. Natl. Acad. Sci. (2025)](https://doi.org/10.1073/pnas.2503439122) | [code](https://github.com/hu-yanxiao/SUS2-MLIP)  

* **Machine Learning Potentials for Alloys: A Detailed Workflow to Predict Phase Diagrams and Benchmark Accuracy** [2025]  
Siya Zhu and Doğuhan Sarıtürk and Raymundo Arroyave.  
[arXiv:2506.16771(2025)](https://arxiv.org/abs/2506.16771) | [code](https://github.com/dogusariturk/PhaseForge)  

* **DistMLIP: A Distributed Inference Platform for Machine Learning Interatomic Potentials** [2025]  
Kevin Han, Bowen Deng, Amir Barati Farimani, Gerbrand Ceder.  
[arXiv:2506.02023 (2025)](https://doi.org/10.48550/arXiv.2506.02023)  

* **Weighted Active Space Protocol for Multireference Machine-Learned Potentials** [2025]  
Aniruddha Seal, Simone Perego, Matthew R. Hennefarth, Umberto Raucci, Luigi Bonati, Andrew L. Ferguson, Michele Parrinello, Laura Gagliardi.  
[arXiv:2505.10505 (2025)](https://arxiv.org/abs/2505.10505)  

* **Discovering High-Entropy Oxides with a Machine-Learning Interatomic Potential** [2025]  
Sivak, Jacob T., Saeed SI Almishal, Mary Kathleen Caucci, Yueze Tan, Dhiya Srikanth, Joseph Petruska, Matthew Furst et al.  
[Phys. Rev. Lett. (2025)](https://doi.org/10.1103/PhysRevLett.134.216101)  

* **Fast and Fourier Features for Transfer Learning of Interatomic Potentials** [2025]  
Pietro Novelli, Giacomo Meanti, Pedro J. Buigues, Lorenzo Rosasco, Michele Parrinello, Massimiliano Pontil, Luigi Bonati.  
[arXiv:2505.05652 (2025)](https://doi.org/10.48550/arXiv.2505.05652) | [code](https://franken.readthedocs.io/)  

* **Practical Machine Learning Strategies. 2. Accurate Prediction of ωB97X-V/6-311+G(2df,2p), ωB97M-V/6-311+G(2df,2p) and ωB97M(2)/6-311+G(2df,2p) Energies From Neural Networks Trained From ωB97X-D/6-31G* Equilibrium Geometries and Energies** [2025]  
Philip Klunzinger, Thomas Hehre, Bernard Deppmeier, William Ohlinger, Warren Hehre.  
[Computer Physics Communications (2025)](https://doi.org/10.1002/jcc.70129)  

* **Structure and Dynamics of CO2 at the Air–Water Interface from Classical and Neural Network Potentials** [2025]  
Nitesh Kumar and Vyacheslav S. Bryantsev.  
[J. Phys. Chem. Lett. (2025)](https://doi.org/10.1021/acs.jpclett.5c01206)  

* **chemtrain: Learning deep potential models via automatic differentiation and statistical physics** [2025]  
Fuchs, Paul, Stephan Thaler, Sebastien Röcken, and Julija Zavadlav.  
[Computer Physics Communications (2025)](https://doi.org/10.1016/j.cpc.2025.109512) | [code](https://github.com/tummfm/chemtrain)  

* **Accurate and efficient machine learning interatomic potentials for finite temperature modelling of molecular crystals** [2025]  
Della Pia, Flaviano and Shi, Benjamin X. and Kapil, Venkat and Zen, Andrea and Alfè, Dario and Michaelides, Angelos.  
[Chem. Sci. (2025)](https://doi.org/10.1039/D5SC01325A) | [code](https://github.com/water-ice-group/MolCrys-MACE)  

* **Impact of Derivative Observations on Gaussian Process Machine Learning Potentials: A Direct Comparison of Three Modeling Approaches** [2025]  
Yulian T. Manchev and Paul L. A. Popelier.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.5c00344)  

* **INN-FF: A Scalable and Efficient Machine Learning Potential for Molecular Dynamics** [2025]  
Taskin Mehereen, Sourav Saha, Intesar Jawad Jaigirdar, Chanwook Park.  
[arXiv:2505.18141 (2025)](https://doi.org/10.48550/arXiv.2505.18141)  

* **Transferability of MACE Graph Neural Network for Range Corrected Δ-Machine Learning Potential QM/MM Applications** [2025]  
Timothy J. Giese, Jinzhe Zeng, and Darrin M. York.  
[J. Phys. Chem. B (2025)](https://doi.org/10.1021/acs.jpcb.5c02006)  

* **OMNI-P2x: A Universal Neural Network Potential for Excited-State Simulations** [2025]  
Martyka M, Tong X-Y, Jankowska J, Dral PO.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-j207x) | [code](https://github.com/dralgroup/omni-p2x)  

* **Interpolation and differentiation of alchemical degrees of freedom in machine learning interatomic potentials** [2025]  
Nam, J., Peng, J. & Gómez-Bombarelli, R.  
[Nat Commun 16, 4350 (2025)](https://doi.org/10.1038/s41467-025-59543-2) | [code](https://github.com/learningmatter-mit/alchemical-mlip)  

* **Uncertainty quantification for neural network potential foundation models** [2025]  
Bilbrey, J.A., Firoz, J.S., Lee, MS. et al..  
[npj Comput Mater 11, 109 (2025)](https://doi.org/10.1038/s41524-025-01572-y) | [code](https://github.com/pnnl/SNAP)  

* **DeePMD-kit v3: A Multiple-Backend Framework for Machine Learning Potentials** [2025]  
Zeng, Jinzhe, Duo Zhang, Anyang Peng, Xiangyu Zhang, Sensen He, Yan Wang, Xinzijian Liu et al.  
[J. Chem. Theory Comput. (2025)](https://arxiv.org/abs/2502.19161) | [arXiv:2502.19161 (2025)](https://arxiv.org/abs/2502.19161) | [code](https://github.com/njzjz/benchmark-dpv3)  

* **Multi-fidelity learning for interatomic potentials: Low-level forces and high-level energies are all you need** [2025]  
Mitchell Messerly, Sakib Matin, Alice E. A. Allen, Benjamin Nebgen, Kipton Barros, Justin S. Smith, Nicholas Lubbers, Richard Messerly.  
[arXiv:2505.01590 (2025)](https://doi.org/10.48550/arXiv.2505.01590)  

* **Egret-1: Pretrained Neural Network Potentials For Efficient and Accurate Bioorganic Simulation** [2025]  
Corin C. Wagen, Elias L. Mann, Jonathon E. Vandezande, Arien M. Wagen, Spencer C. Schneider.  
[arXiv:2504.20955(2025)](https://arxiv.org/abs/2504.20955) | [code](https://github.com/rowansci/egret-public)  

* **AIMNet2: A Neural Network Potential to Meet your Neutral, Charged, Organic, and Elemental-Organic Needs** [2025]  
Isayev, Olexandr and Anstine, Dylan and Zubaiuk, Roman.  
[Chem. Sci. (2025)](http://dx.doi.org/10.1039/D4SC08572H) | [code](https://github.com/isayevlab/aimnetcentral)  

* **Advancing Density Functional Tight-Binding method for Large Organic Molecules through Equivariant Neural Networks** [2025]  
Medrano Sandonas LR, Puleva M, Parra Payano R, Stöhr M, Cuniberti G, Tkatchenko A.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-z3mhh) | [code](https://github.com/lmedranos/EquiDTB)  

* **Evaluating Universal Interatomic Potentials for Molecular Dynamics of Real-World Minerals** [2025]  
Sajid Mannan, Carmelo Gonzales, Vaibhav Bihani, Kin Long Kelvin Lee, Nitya Nand Gosvami, Santiago Miret, N M Anoop Krishnan.  
[ ICLR 2025 Workshop AI4MAT (2025)](https://openreview.net/forum?id=me0flBb1hi)  

* **Improving robustness and training efficiency of machine-learned potentials by incorporating short-range empirical potentials** [2025]  
Yan, Zihan, Zheyong Fan, and Yizhou Zhu.  
[arXiv:2504.15925 (2025)](https://doi.org/10.48550/arXiv.2504.15925) | [code](https://github.com/zhyan0603/GPUMDkit)  

* **Automated Pynta-Based Curriculum for ML-Accelerated Calculation of Transition States** [2025]  
revor Price, Saurabh Sivakumar, Matthew S. Johnson, Judit Zádor, and Ambarish Kulkarni.  
[J. Phys. Chem. C (2025)](https://doi.org/10.1021/acs.jpcc.5c00305) | [code](https://github.com/tdprice-858/Pynta-ML)  

* **Extending atomic decomposition and many-body representation with a chemistry-motivated approach to machine learning potentials** [2025]  
Yu, Q., Ma, R., Qu, C. et al.  
[Nat Comput Sci (2025)](https://doi.org/10.1038/s43588-025-00790-0) | [code](https://doi.org/10.5281/zenodo.14954863)  

* **Learning Pairwise Interaction for Extrapolative and Interpretable Machine Learning Interatomic Potentials with Physics-Informed Neural Network** [2025]  
Hoje Chun, Minjoon Hong, Seung Hyo Noh, and Byungchan Han.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.5c00090) | [code](https://figshare.com/s/d13b57f2fbca1d0d4b57)  

* **Accurate and Affordable Simulation of Molecular Infrared Spectra with AIQM Models** [2025]  
Yi-Fan Hou, Cheng Wang, and Pavlo O. Dral.  
[J. Phys. Chem. A  (2025)](https://doi.org/10.1021/acs.jpca.5c00146) | [code](https://github.com/dralgroup/mlatom)  

* **PET-MAD, a universal interatomic potential for advanced materials modeling** [2025]  
BArslan Mazitov, Filippo Bigi, Matthias Kellner, Paolo Pegolo, Davide Tisi, Guillaume Fraux, Sergey Pozdnyakov, Philip Loche, Michele Ceriotti.  
[arXiv:2503.14118 (2025)](https://doi.org/10.48550/arXiv.2503.14118) | [code](https://github.com/lab-cosmo/pet-mad)  

* **Generator of Neural Network Potential for Molecular Dynamics: Constructing Robust and Accurate Potentials with Active Learning for Nanosecond-Scale Simulations** [2025]  
Benjamin Rhodes, Sander Vandenhaute, Vaidotas Šimkus, James Gin, Jonathan Godwin, Tim Duignan, Mark Neumann.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.4c01613)  

* **Orb-v3: atomistic simulation at scale** [2025]  
Benjamin Rhodes, Sander Vandenhaute, Vaidotas Šimkus, James Gin, Jonathan Godwin, Tim Duignan, Mark Neumann.  
[arXiv:2504.06231 (2025)](https://arxiv.org/abs/2504.06231) | [code](https://github.com/orbital-materials/orb-models)  

* **Orb-v3: atomistic simulation at scale** [2025]  
Benjamin Rhodes, Sander Vandenhaute, Vaidotas Šimkus, James Gin, Jonathan Godwin, Tim Duignan, Mark Neumann.  
[arXiv:2504.06231 (2025)](https://arxiv.org/abs/2504.06231) | [code](https://github.com/orbital-materials/orb-models)  

* **Carbonic anhydrase II simulated with a universal neural network potential** [2025]  
Timothy T. Duignan.  
[arXiv:2503.13789 (2025)](https://arxiv.org/abs/2503.13789) | [code](https://github.com/timduignan/CAII-Orb-d3-v2)  

* **Machine learning interatomic potential can infer electrical response** [2025]  
Peichen Zhong, Dongjin Kim, Daniel S. King, Bingqing Cheng.  
[arXiv:2504.05169 (2025)](https://arxiv.org/abs/2504.05169) | [code](https://github.com/BingqingCheng/cace)  

* **QuantumBind-RBFE: Accurate Relative Binding Free Energy Calculations Using Neural Network Potentials** [2025]  
Zariquiey, Francesc Sabanés, Stephen E. Farr, Stefan Doerr, and Gianni De Fabritiis.  
[J. Chem. Inf. Model. (2025)](https://doi.org/10.1021/acs.jcim.5c00033) | [arXiv:2501.01811 (2025)](https://arxiv.org/abs/2501.01811) | [code](https://github.com/Acellera/quantumbind_rbfe)  

* **DeePMD-GNN: A DeePMD-kit Plugin for External Graph Neural Network Potentials** [2025]  
Jinzhe Zeng, Timothy J. Giese, Duo Zhang, Han Wang, and Darrin M. York.  
[J. Chem. Inf. Model. (2025)](https://doi.org/10.1021/acs.jcim.4c02441) | [code](https://gitlab.com/RutgersLBSR/deepmd-gnn)  

* **An Investigation of Physics Informed Neural Networks to Solve the Poisson–Boltzmann Equation in Molecular Electrostatics** [2025]  
Martín A. Achondo, Jehanzeb H. Chaudhry, and Christopher D. Cooper.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.4c01747) | [code](https://github.com/MartinAchondo/XPPBE)  

* **Hierarchical Deep Potential with Structure Constraints for Efficient Coarse-Grained Modeling** [2025]  
Qi Huang, Yedi Li, Lei Zhu, and Wenjie Yu.  
[J. Chem. Inf. Model.(2025)](https://doi.org/10.1021/acs.jcim.4c02042) | [code](https://github.com/huang-qi/HDP-SC/)  

* **Does Hessian Data Improve the Performance of Machine Learning Potentials?** [2025]  
Austin Rodriguez, Justin S. Smith, Jose L. Mendoza-Cortes.  
[arXiv:2503.07839 (2025)](https://arxiv.org/abs/2503.07839)  

* **Free energy profiles for chemical reactions in solution from high-dimensional neural network potentials: The case of the Strecker synthesis** [2025]  
Alea Miako Tokita, Timothée Devergne, A. Marco Saitta, Jörg Behler.  
[arXiv:2503.05370 (2025)](https://arxiv.org/abs/2503.05370)  

* **Accurate Free Energy Calculation via Multiscale Simulations Driven by Hybrid Machine Learning and Molecular Mechanics Potentials** [2025]  
Wang X, Wu X, Brooks B, Wang J.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2024-zq975-v2)  

* **ANI-1xBB: an ANI based reactive potential** [2025]  
Zhang S, Zubatyuk R, Yang Y, Roitberg A, Isayev O.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-m2nqq)  

* **Efficient Training of Neural Network Potentials for Chemical and Enzymatic Reactions by Continual Learning** [2025]  
Yao-Kun Lei, Kiyoshi Yagi, and Yuji Sugita.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.4c01393)   

* **Efficient equivariant model for machine learning interatomic potentials** [2025]  
Yang, Z., Wang, X., Li, Y. et al.  
[npj Comput Mater 11, 49 (2025)](https://doi.org/10.1038/s41524-025-01535-3) | [code](https://github.com/Shen-Group/E2GNN)   

* **Neural Network Potential with Multiresolution Approach Enables Accurate Prediction of Reaction Free Energies in Solution** [2025]  
Felix Pultar, Moritz Thürlemann, Igor Gordiy, Eva Doloszeski, and Sereina Riniker.  
[J. Am. Chem. Soc. (2025)](https://doi.org/10.1021/jacs.4c17015) | [code](https://github.com/rinikerlab/amp_qmmm)  

* **End-To-End Learning of Classical Interatomic Potentials for Benchmarking Anion Polarization Effects in Lithium Polymer Electrolytes** [2025]  
Pablo A. Leon, Avni Singhal, Jurgis Ruza, Jeremiah A. Johnson, Yang Shao-Horn, and Rafael Gomez-Bombarelli.  
[Chem. Mater. (2025)](https://doi.org/10.1021/acs.chemmater.4c02529) | [code](https://github.com/learningmatter-mit/AutoBADDIE)   

* **Learning Smooth and Expressive Interatomic Potentials for Physical Property Prediction** [2025]  
Xiang Fu, Brandon M. Wood, Luis Barroso-Luque, Daniel S. Levine, Meng Gao, Misko Dzamba, C. Lawrence Zitnick.  
[arXiv:2502.12147 (2025)](https://arxiv.org/abs/2502.12147) | [code](https://github.com/FAIR-Chem/fairchem)   

* **PiNN: Equivariant Neural Network Suite for Modeling Electrochemical Systems** [2025]  
Jichen Li, Lisanne Knijff, Zhan-Yun Zhang, Linnéa Andersson, and Chao Zhang.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.4c01570) | [code](https://teoroo-cmc.github.io/PiNN/master/usage/pinet)   

* **An automated framework for exploring and learning potential-energy surfaces** [2024]  
Liu, Yuanbin, Joe D. Morrow, Christina Ertural, Natascia L. Fragapane, John LA Gardner, Aakash A. Naik, Yuxing Zhou, Janine George, and Volker L. Deringer.  
[arXiv:2412.16736 (2024)](https://arxiv.org/abs/2412.16736) | [code](https://github.com/autoatml/autoplex)   

* **The Importance of Being Scalable: Improving the Speed and Accuracy of Neural Network Interatomic Potentials Across Chemical Domains** [2024]  
Eric Qu, Aditi S. Krishnapriyan.  
[arXiv:2410.24169 (2024)](https://arxiv.org/abs/2410.24169v1)   

* **Cartesian atomic moment machine learning interatomic potentials** [2024]  
Mingjian Wen, Wei-Fan Huang, Jin Dai, Santosh Adhikari.  
[arXiv:2411.12096 (2024)](https://doi.org/10.48550/arXiv.2411.12096) | [code](https://github.com/wengroup/camp)  

* **DeepConf: Leveraging ANI-ML Potentials for Exploring Local Minima with A Focus on Bioactive Conformations** [2024]  
Tayfuroglu O, Zengin İN, Koca MS, Kocak A.  
[ChemRxiv. (2024)](https://doi.org/10.26434/chemrxiv-2024-xqbdp-v2) | [code](https://github.com/otayfuroglu/DeepConf)  

* **Universal Machine Learning Interatomic Potentials are Ready for Phonons** [2024]  
Antoine Loew, Dewen Sun, Hai-Chen Wang, Silvana Botti, Miguel A. L. Marques.  
[arXiv:2412.16551 (2024)](https://arxiv.org/abs/2412.16551) | [code](https://github.com/hyllios/utils)   

* **Neural Network Potentials for Enabling Advanced Small-Molecule Drug Discovery and Generative Design** [2024]  
Barnett, Simon, and John D. Chodera.  
[GEN Biotechnology 3.3 (2024)](https://doi.org/10.1089/genbio.2024.0011)    

* **Online Test-time Adaptation for Interatomic Potentials** [2024]  
Taoyong Cui, Chenyu Tang, Dongzhan Zhou, Yuqiang Li, Xingao Gong, Wanli Ouyang, Mao Su, Shufei Zhang.  
[arXiv:2405.08308 (2024)](https://arxiv.org/abs/2405.08308) | [code](https://github.com/TaoyongCui/TAIP-codes)   

* **Uncertainty-biased molecular dynamics for learning uniformly accurate interatomic potentials** [2024]  
Zaverkin, V., Holzmüller, D., Christiansen, H. et al.  
[npj Comput Mater 10, 83 (2024)](https://doi.org/10.1038/s41524-024-01254-1) | [code](https://github.com/nec-research/alebrew)   

* **General-purpose machine-learned potential for 16 elemental metals and their alloys** [2024]  
Song, K., Zhao, R., Liu, J. et al.  
[Nat Commun 15, 10208 (2024)](https://doi.org/10.1038/s41467-024-54554-x) | [code](https://github.com/brucefan1983/GPUMD)    

* **Scalable Parallel Algorithm for Graph Neural Network Interatomic Potentials in Molecular Dynamics Simulations** [2024]  
Yutack Park, Jaesun Kim, Seungwoo Hwang, and Seungwu Han.  
[J. Chem. Theory Comput. (2024)](https://doi.org/10.1021/acs.jctc.4c00190) | [code](https://github.com/MDIL-SNU/SevenNet)  

* **Ab initio Accuracy Neural Network Potential for Drug-like Molecules** [2024]  
Yang M, Zhang D, Wang X, Zhang L, Zhu T, Wang H.  
[ChemRxiv. (2024)](https://doi.org/10.26434/chemrxiv-2024-sq8nh) | [data](https://www.aissquare.com/datasets/detail?pageType=datasets&name=Drug%28drug-like-molecule%29_DPA_v1_0&id=143)    

* **HH130: a standardized database of machine learning interatomic potentials, datasets, and its applications in the thermal transport of half-Heusler thermoelectrics** [2024]  
Yang, Yuyan, Yifei Lin, Shengnan Dai, Yifan Zhu, Jinyang Xi, Lili Xi, Xiaokun Gu, David J. Singh, Wenqing Zhang, and Jiong Yang.  
[Digital Discovery (2024)](https://doi.org/10.1039/D4DD00240G) | [data](http://www.mathub3d.net)    

* **Efficient Training of Neural Network Potentials for Chemical and Enzymatic Reactions by Continual Learning** [2024]  
Lei Y-K, Yagi K, Sugita Y.  
[ChemRxiv. (2024)](https://doi.org/10.26434/chemrxiv-2024-xkxd5)  

* **Learning together: Towards foundation models for machine learning interatomic potentials with meta-learning** [2024]  
Allen, A.E.A., Lubbers, N., Matin, S. et al.  
[npj Comput Mater 10, 154 (2024)](https://doi.org/10.1038/s41524-024-01339-x)  

* **Quantum-accurate machine learning potentials for metal-organic frameworks using temperature driven active learning** [2024]  
Sharma, A., Sanvito, S.  
[npj Comput Mater 10, 237 (2024)](https://doi.org/10.1038/s41524-024-01427-y) | [code](https://github.com/asharma-ms/MOF_MLP_2024)  

* **Molecular Simulations with a Pretrained Neural Network and Universal Pairwise Force Fields** [2024]  
Kabylda A, Frank JT, Dou SS, Khabibrakhmanov A, Sandonas LM, Unke OT, et al.  
[ChemRxiv. (2024)](https://doi.org/10.26434/chemrxiv-2024-bdfr0) | [code](https://github.com/general-molecular-simulations/so3lr)  

* **AMARO: All Heavy-Atom Transferable Neural Network Potentials of Protein Thermodynamics** [2024]  
Mirarchi, Antonio, Raul P. Pelaez, Guillem Simeon, and Gianni De Fabritiis.  
[arXiv:2409.17852 (2024)](https://arxiv.org/abs/2409.17852) | [code](https://github.com/compsciencelab/amaro)  

* **Revisiting Aspirin Polymorphic Stability Using a Machine Learning Potential** [2024]  
Hattori, Shinnosuke, and Qiang Zhu.  
[ACS Omega (2024)](https://doi.org/10.1021/acsomega.4c04782)  

* **Enhancing Protein–Ligand Binding Affinity Predictions Using Neural Network Potentials** [2024]  
Sabanés Zariquiey, F., Galvelis, R., Gallicchio, E., Chodera, J.D., Markland, T.E. and De Fabritiis, G.   
[J. Chem. Inf. Model. (2024)](https://doi.org/10.1021/acs.jcim.3c02031) | [code](https://github.com/compsciencelab/ATM_benchmark/tree/main/ATM_With_NNPs)  

* **Universal-neural-network-potential molecular dynamics for lithium metal and garnet-type solid electrolyte interface** [2024]  
Iwasaki, R., Tanibata, N., Takeda, H. et al.   
[Commun Mater 5, 148 (2024)](https://doi.org/10.1038/s43246-024-00595-0)  

* **The Potential of Neural Network Potentials** [2024]  
Duignan, Timothy T.   
[ACS Physical Chemistry Au 4.3 (2024)](https://doi.org/10.1021/acsphyschemau.4c00004)  

* **GPIP: Geometry-enhanced Pre-training on Interatomic Potentials** [2024]  
Su, M., S. Zhang, T. Cui, C. Tang, Y. Li, Y. Dong, X. Gong, W. Ouyang, and L. Bai.  
[arXiv:2309.15718 (2024)](https://doi.org/10.1038/s42256-024-00818-6) | [code](https://github.com/cuitaoyong/GPIP)  

* **AIMNet2: A Neural Network Potential to Meet your Neutral, Charged, Organic, and Elemental-Organic Needs** [2024]  
Anstine, Dylan, Roman Zubatyuk, and Olexandr Isayev.   
[chemrxiv-2023-296ch-v2 (2024)](https://doi.org/10.26434/chemrxiv-2023-296ch-v2) | [code](https://github.com/isayevlab/aimnet2)  

* **NNP/MM: Accelerating Molecular Dynamics Simulations with Machine Learning Potentials and Molecular Mechanics** [2023]  
Galvelis, R., Varela-Rial, A., Doerr, S., Fino, R., Eastman, P., Markland, T.E., Chodera, J.D. and De Fabritiis, G.   
[J. Chem. Inf. Model. (2023)](https://doi.org/10.1021/acs.jcim.3c00773) | [code](https://github.com/openmm/nnpops)  

* **CHGNet as a pretrained universal neural network potential for charge-informed atomistic modelling** [2023]  
Deng, B., Zhong, P., Jun, K. et al.  
[Nat Mach Intell 5, 1031–1041 (2023)](https://doi.org/10.1038/s42256-023-00716-3) | [code](https://github.com/CederGroupHub/chgnet)  

* **Towards universal neural network potential for material discovery applicable to arbitrary combination of 45 elements** [2022]  
Takamoto, S., Shinagawa, C., Motoki, D. et al.   
[Nat Commun 13, 2991 (2022)](https://doi.org/10.1038/s41467-022-30687-9) | [data](https://doi.org/10.6084/m9.figshare.19658538)  

* **E(3)-equivariant graph neural networks for data-efficient and accurate interatomic potentials** [2022]  
Batzner, S., Musaelian, A., Sun, L. et al.   
[Nat Commun 13, 2453 (2022)](https://doi.org/10.1038/s41467-022-29939-5) | [data](https://github.com/teddykoker/nequip-eqx)  

* **Teaching a neural network to attach and detach electrons from molecules** [2021]  
Zubatyuk, R., Smith, J.S., Nebgen, B.T. et al.   
[Nat Commun 12, 4870 (2021)](https://doi.org/10.1038/s41467-021-24904-0)  | [code](https://github.com/isayevlab/aimnetnse)  

* **Four Generations of High-Dimensional Neural Network Potentials** [2021]  
Behler, Jorg.   
[Chemical Reviews 121.16 (2021)](https://doi.org/10.1021/acs.chemrev.0c00868)  

* **DP-GEN: A concurrent learning platform for the generation of reliable deep learning based potential energy models** [2020]  
Zhang, Yuzhi, Haidi Wang, Weijie Chen, Jinzhe Zeng, Linfeng Zhang, Han Wang, and E. Weinan.   
[Computer Physics Communications 253 (2020)](https://doi.org/10.1016/j.cpc.2020.107206)  | [code](https://github.com/deepmodeling/dpgen)  










### Neural Reactive Potential




* **AIMNet2-NSE: A Transferable Reactive Neural Network Potential for Open-Shell Chemistry** [2025]  
Kalita B, Zubatyuk R, Anstine DM, Bergeler M, Settels V, Stork C, et al.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-kdg6n) | [data](https://figshare.com/s/988db758cceb66f220c1?file=56714804)  

* **Machine learning prediction of a chemical reaction over 8 decades of energy** [2025]  
Daniel Julian, Jesús Pérez-Ríos.  
[arXiv:2507.01793 (2025)](https://doi.org/10.48550/arXiv.2507.01793)  

* **AIMNet2-rxn: A Machine Learned Potential for Generalized Reaction Modeling on a Millions-of-Pathways Scale** [2025]  
Anstine DM, Zhao Q, Zubatiuk R, Zhang S, Singla V, Nikitin F, et al.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-hpdmg) | [code](https://github.com/isayevlab/aimnetcentral)  

* **Harnessing Machine Learning to Enhance Transition State Search with Interatomic Potentials and Generative Models** [2025]   
 Zhao Q, Han Y, Zhang D, Wang J, Zhong P, Cui T, et al.  
  [ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-mt6hc-v2)  

* **Capturing Excited State Proton Transfer Dynamics with Reactive Machine Learning Potentials** [2025]   
 Umberto Raucci.  
  [J. Phys. Chem. Lett. (2025)](https://doi.org/10.1021/acs.jpclett.5c00688)  

* **ANI-1xBB: An ANI-Based Reactive Potential for Small Organic Molecules** [2025]   
 Shuhao Zhang, Roman Zubatyuk, Yinuo Yang, Adrian Roitberg, and Olexandr Isayev.  
  [J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.5c00347) | [code](https://github.com/amateurcat/ANI-1xBB)  

* **The evolution of machine learning potentials for molecules, reactions and materials** [2025]   
 Xia, Junfan and Zhang, Yaolong and Jiang, Bin.  
  [Chem. Soc. Rev. (2025)](https://doi.org/10.1039/D5CS00104H)  

* **Reaction dynamics of Diels–Alder reactions from machine learned potentials** [2022]  
 Young, Tom A., Tristan Johnston-Wood, Hanwen Zhang, and Fernanda Duarte.  
  [Physical Chemistry Chemical Physics 24.35 (2022)](https://doi.org/10.1039/D2CP02978B) | [code](https://github.com/duartegroup/mlp-train)  








### Reactive Force Fields





* **Boosting ReaxFF Reactive Force Field Optimization with Adaptive Sampling** [2025]   
 Shuang Li, Siyuan Yang, Sibing Chen, Wei Zheng, Zejian Dong, Langli Luo, Weiwei Zhang, and Xing Chen.  
  [J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.4c01748) | [code](https://github.com/tju-chen-group/Adaptive-Sampling-Optimization-code)  








### Free Energy Perturbation





* **Machine Learning Guided AQFEP: A Fast and Efficient Absolute Free Energy Perturbation Solution for Virtual Screening** [2025]  
Jordan E. Crivelli-Decker, Zane Beckwith, Gary Tom, Ly Le, Sheenam Khuttan, Romelia Salomon-Ferrer, Jackson Beall, Rafael Gómez-Bombarelli, and Andrea Bortolato.  
  [J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.4c00399) | [data](https://zenodo.org/records/12827945)  

* **Considerations in the use of machine learning force fields for free energy calculations** [2025]  
Orlando A. Mendible-Barreto, Jonathan K. Whitmer, Yamil J. Colón.  
  [J. Chem. Phys. (2025)](https://doi.org/10.1063/5.0252043) | [code](https://github.com/omendibleba/Considerations_for_MLPs_FES)  

* **Narjes Ansari, Zhifeng Francis Jing, Antoine Gagelin, Florent Hédin, Félix Aviat, Jérôme Hénin, Jean-Philip Piquemal, and Louis Lagardère** [2025]  
Narjes Ansari, Zhifeng Francis Jing, Antoine Gagelin, Florent Hédin, Félix Aviat, Jérôme Hénin, Jean-Philip Piquemal, and Louis Lagardère.  
  [J. Phys. Chem. Lett. (2025)](https://doi.org/10.1021/acs.jpclett.5c00683) | [code](https://github.com/ansarinarjes/ABFE-L-ABF-OPES)  

* **Applying Absolute Free Energy Perturbation Molecular Dynamics to Diffusively Binding Ligands** [2025]  
Laracuente, Xavier E., Bryan M. Delfing, Xingyu Luo, Audrey Olson, William Jeffries, Steven R. Bowers, Kenneth W. Foreman et al.  
  [J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.5c00121) | [code](https://github.com/KlimovLab/AFEP-REST_KKPK_impa)  

* **QuantumBind-RBFE: Accurate Relative Binding Free Energy Calculations Using Neural Network Potentials** [2025]  
Zariquiey, Francesc Sabanés, Stephen E. Farr, Stefan Doerr, and Gianni De Fabritiis.  
[J. Chem. Inf. Model. (2025)](https://doi.org/10.1021/acs.jcim.5c00033) | [arXiv:2501.01811 (2025)](https://arxiv.org/abs/2501.01811) | [code](https://github.com/Acellera/quantumbind_rbfe)  

* **Free energy profiles for chemical reactions in solution from high-dimensional neural network potentials: The case of the Strecker synthesis** [2025]  
Alea Miako Tokita, Timothée Devergne, A. Marco Saitta, Jörg Behler.  
[arXiv:2503.05370 (2025)](https://arxiv.org/abs/2503.05370)  

* **Accurate Free Energy Calculation via Multiscale Simulations Driven by Hybrid Machine Learning and Molecular Mechanics Potentials** [2025]  
Wang X, Wu X, Brooks B, Wang J.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2024-zq975-v2)  

* **Robust protein–ligand interaction modeling through integrating physical laws and geometric knowledge for absolute binding free energy calculation** [2025]   
Su, Qun, Jike Wang, Qiaolin Gou, Renling Hu, Linlong Jiang, Hui Zhang, Tianyue Wang et al.  
  [Chemical Science (2025)](https://doi.org/10.1039/D4SC07405J) | [code](https://github.com/lingcon01/LumiNet)  

* **Lambda-ABF-OPES: Faster Convergence with High Accuracy in Alchemical Free Energy Calculations** [2025]  
Narjes Ansari, Francis Jing, Antoine Gagelin, Florent Hédin, Félix Aviat, Jérôme Hénin, Jean-Philip Piquemal, Louis Lagardère.  
  [arXiv:2502.17233 (2025)](https://arxiv.org/abs/2502.17233)   

* **Comparison of Methodologies for Absolute Binding Free Energy Calculations of Ligands to Intrinsically Disordered Proteins** [2024]   
Michail Papadourakis, Zoe Cournia, Antonia S. J. S. Mey, and Julien Michel.    
  [J. Chem. Theory Comput. (2024)](https://doi.org/10.1021/acs.jctc.4c00942) | [code](https://github.com/michellab/idpabfe)  

* **FEP-SPell-ABFE: An Open-Source Automated Alchemical Absolute Binding Free Energy Calculation Workflow for Drug Discovery** [2024]   
Pengfei Li,Tingting Pu ,Ye Mei.    
  [ChemRxiv. (2024)](https://doi.org/10.26434/chemrxiv-2024-tkvrh) | [code](https://github.com/freeenergylab/FEP-SPell-ABFE)    

* **Studying the Collective Functional Response of a Receptor in Alchemical Ligand Binding Free Energy Simulations with Accelerated Solvation Layer Dynamics** [2024]   
Wei Jiang.   
  [J. Chem. Theory Comput. (2024)](https://doi.org/10.1021/acs.jctc.4c00191)  

* **Alchemical Transformations and Beyond: Recent Advances and Real-World Applications of Free Energy Calculations in Drug Discovery** [2024]   
 Qian, Runtong, Jing Xue, You Xu, and Jing Huang.   
  [J. Chem. Inf. Model. (2024)](https://doi.org/10.1021/acs.jcim.4c01024)  

* **Automated Adaptive Absolute Binding Free Energy Calculations** [2024]  
Clark, Finlay, Graeme Robb, Daniel Cole, and Julien Michel.   
[J. Chem. Theory Comput. (2024)](https://doi.org/10.1021/acs.jctc.4c00806) | [code](https://github.com/michellab/a3fe)  

* **Machine Learning Guided AQFEP: A Fast and Efficient Absolute Free Energy Perturbation Solution for Virtual Screening** [2024]  
Crivelli-Decker, J.E., Beckwith, Z., Tom, G., Le, L., Khuttan, S., Salomon-Ferrer, R., Beall, J., Gómez-Bombarelli, R. and Bortolato, A.   
[J. Chem. Theory Comput. (2024)](https://doi.org/10.1021/acs.jctc.4c00399) | [code](https://zenodo.org/records/12827945)  

* **The maximal and current accuracy of rigorous protein-ligand binding free energy calculations** [2023]   
Ross, G.A., Lu, C., Scarabelli, G. et al.    
  [Commun Chem 6, 222 (2023)](https://doi.org/10.1038/s42004-023-01019-9) | [code](https://github.com/schrodinger/public_binding_free_energy_benchmark)    







### Ab Initio










* **LAMBench: A Benchmark for Large Atomistic Models** [2025]  
Anyang Peng, Chun Cai, Mingyu Guo, Duo Zhang, Chengqian Zhang, Wanrun Jiang, Yinan Wang, Antoine Loew, Chengkun Wu, Weinan E, Linfeng Zhang, Han Wang.  
[arXiv:2504.19578 (2025)](https://doi.org/10.48550/arXiv.2504.19578) | [code](https://github.com/deepmodeling/lambench)  

* **DREAMS: Density Functional Theory Based Research Engine for Agentic Materials Simulation** [2025]  
Ziqi Wang, Hongshuo Huang, Hancheng Zhao, Changwen Xu, Shang Zhu, Jan Janssen, Venkatasubramanian Viswanathan.  
[arXiv:2507.14267 (2025)](https://doi.org/10.48550/arXiv.2507.14267) | [code](https://github.com/BattModels/material_agent)  

* **Transferring Knowledge from MM to QM: A Graph Neural Network-Based Implicit Solvent Model for Small Organic Molecules** [2025]  
Xu, M., Wang, S., He, Y. et al.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.5c00728) | [code](https://github.com/rinikerlab/QM-GNNIS)  

* **Efficient modeling of ionic and electronic interactions by a resistive memory-based reservoir graph neural network** [2025]  
Xu, M., Wang, S., He, Y. et al.  
[Nat Comput Sci (2025)](https://doi.org/10.1038/s43588-025-00844-3) | [code](https://github.com/hustmeng/RGNN)  

* **OrbitAll: A Unified Quantum Mechanical Representation Deep Learning Framework for All Molecular Systems** [2025]  
Kang, Beom Seok, Vignesh C. Bhethanabotla, Amin Tavakoli, Maurice D. Hanisch, William A. Goddard III, and Anima Anandkumar.  
[arXiv:2507.03853 (2025)](https://doi.org/10.48550/arXiv.2507.03853)  

* **An ab initio foundation model of wavefunctions that accurately describes chemical bond breaking** [2025]  
Adam Foster, Zeno Schätzle, P. Bernát Szabó, Lixue Cheng, Jonas Köhler, Gino Cassella, Nicholas Gao, Jiawei Li, Frank Noé, Jan Hermann.  
[arXiv:2506.19960 (2025)](https://doi.org/10.48550/arXiv.2506.19960) | [code](https://github.com/microsoft/oneqmc)  

* **Predicting Oxidation Potentials with DFT-Driven Machine Learning** [2025]  
Shweta Sharma, Natan Kaminsky, Kira Radinsky, and Lilac Amirav.  
[J. Chem. Inf. Model. (2025)](https://doi.org/10.1021/acs.jcim.5c00159) | [code](https://github.com/nkami/oxpot)  

* **g-xTB: A General-Purpose Extended Tight-Binding Electronic Structure Method For the Elements H to Lr (Z=1–103)** [2025]  
Froitzheim T, Müller M, Hansen A, Grimme S.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-bjxvt) | [code](https://github.com/grimme-lab/g-xtb)  

* **Discovery of chemically modified higher tungsten boride by means of hybrid GNN/DFT approach** [2025]  
Matsokin, N.A., Eremin, R.A., Kuznetsova, A.A. et al.  
[npj Comput Mater 11, 163 (2025)](https://doi.org/10.1038/s41524-025-01628-z)  

* **Revisiting a large and diverse data set for barrier heights and reaction energies: best practices in density functional theory calculations for chemical kinetics** [2025]  
Liu, Xiao and Spiekermann, Kevin A. and Menon, Angiras and Green, William H. and Head-Gordon, Martin.  
[Phys. Chem. Chem. Phys. (2025)](https://doi.org/10.1039/D5CP01181G)  

* **Accurate and scalable exchange-correlation with deep learning** [2025]  
Luise, Giulia et al.  
[arXiv:2506.14665 (2025)](https://doi.org/10.48550/arXiv.2506.14665) | [data](https://zenodo.org/records/15387280)  

* **Self-Refining Training for Amortized Density Functional Theory** [2025]  
Majdi Hassan, Cristian Gabellini, Hatem Helal, Dominique Beaini, Kirill Neklyudov.  
[arXiv:2506.01225 (2025)](https://arxiv.org/abs/2506.01225) | [code](https://github.com/majhas/self-refining-dft)  

* **Unified deep learning framework for many-body quantum chemistry via Green’s functions** [2025]  
SVenturella, C., Li, J., Hillenbrand, C. et al.  
[Nat Comput Sci (2025)](https://doi.org/10.1038/s43588-025-00810-z) | [code](https://github.com/ZhuGroup-Yale/mlgf)  

* **High-order Equivariant Flow Matching for Density Functional Theory Hamiltonian Prediction** [2025]  
Seongsu Kim, Nayoung Kim, Dongwoo Kim, Sungsoo Ahn.  
[arXiv:2505.13424 (2025)](https://doi.org/10.48550/arXiv.2505.13424)  

* **The Enduring Relevance of Semiempirical Quantum Mechanics** [2025]  
Jonathan E. Moussa.  
[arXiv:2505.18817 (2025)](https://doi.org/10.48550/arXiv.2505.18817)  

* **DENSE SENSE : A novel approach utilizing an electron density augmented machine learning paradigm to understand a complex odour landscape** [2025]  
Saha P, Sharma M, Balaji S, Barsainyan AA, Kumar R, Steuber V, et al.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-68d6w-v3) | [code](https://github.com/CSIO-FPIL/Dense-Sense-ODOR)  

* **Advancing Density Functional Tight-Binding method for Large Organic Molecules through Equivariant Neural Networks** [2025]  
Medrano Sandonas LR, Puleva M, Parra Payano R, Stöhr M, Cuniberti G, Tkatchenko A.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-z3mhh) | [code](https://github.com/lmedranos/EquiDTB)  

* **Accurate Electrostatics for Biomolecular Systems through Machine Learning** [2025]  
Hosseini AN, Kriz K, van der Spoel D.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-q2v9m) | [code](https://github.com/AlexandriaChemistry/ACT)  

* **Accurate Ab-initio Neural-network Solutions to Large-Scale Electronic Structure Problems** [2025]  
Michael Scherbela, Nicholas Gao, Philipp Grohs, Stephan Günnemann.  
[arXiv:2504.06087(2025)](https://arxiv.org/abs/2504.06087) | [code](https://arxiv.org/abs/2504.06087)  

* **Analytical ab initio hessian from a deep learning potential for transition state optimization** [2024]  
KYuan, E.CY., Kumar, A., Guan, X. et al.  
[Nat Commun 15, 8865 (2024)](https://doi.org/10.1038/s41467-024-52481-5) | [code](https://github.com/zadorlab/sella)  








## Deep Learning-molecular conformations





* **Accurate Electrostatics for Biomolecular Systems through Machine Learning** [2025]  
Hosseini AN, Kriz K, van der Spoel D.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-q2v9m) | [code](https://github.com/AlexandriaChemistry/ACT)  

* **Molecular Simulations with a Pretrained Neural Network and Universal Pairwise Force Fields** [2024]  
Kabylda A, Frank JT, Dou SS, Khabibrakhmanov A, Sandonas LM, Unke OT, et al.  
[ChemRxiv. (2024)](https://doi.org/10.26434/chemrxiv-2024-bdfr0) | [code](https://github.com/general-molecular-simulations/so3lr)  

* **SpaiNN: Equivariant Message Passing for Excited-State Nonadiabatic Molecular Dynamics** [2024]  
Mausenberger, Sascha, Carolin Müller, Alexandre Tkatchenko, Philipp Marquetand, Leticia González, and Julia Westermayr.   
[Chemical Science (2024)](https://doi.org/10.1039/D4SC04164J) | [code](https://github.com/schnarc/SchNarc)  

* **GLOW: A Workflow Integrating Gaussian-Accelerated Molecular Dynamics and Deep Learning for Free Energy Profiling** [2022]  
Do, Hung N., Jinan Wang, Apurba Bhattarai, and Yinglong Miao.   
[J. Chem. Theory Comput. (2022)](https://doi.org/10.1021/acs.jctc.1c01055) | [code](https://github.com/MiaoLab20/GLOW)  










### AlphaFold-based




* **Applied causality to infer protein dynamics and kinetics** [2025]  
Akashnathan Aranganathan, Eric R. Beyerle.  
[arXiv:2508.12060 (2025)](https://doi.org/10.48550/arXiv.2508.12060)  

* **af2rave: Protein Ensemble Generation with Physics-Based Sampling** [2025]  
Teng D, Meraz VJ, Aranganathan A, Gu X, Tiwary P.  
[Digital Discovery (2025)](https://doi.org/10.1039/D5DD00201J) | [ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-q3mwr) | [code](https://github.com/tiwarylab/af2rave)  

* **Modeling protein conformational ensembles by guiding AlphaFold2 with Double Electron Electron Resonance (DEER) distance distributions** [2025]  
Wu, T., Stein, R.A., Kao, TY. et al.  
[Nat Commun 16, 7107 (2025)](https://doi.org/10.1038/s41467-025-62582-4) | [code](https://github.com/CAAIPD/DEERFold)  

* **Large-scale predictions of alternative protein conformations by AlphaFold2-based sequence association** [2025]  
Lee, M., Schafer, J.W., Prabakaran, J. et al.  
[Nat Commun 16, 5622 (2025)](https://doi.org/10.1021/acs.jcim.5c00489) | [code](https://github.com/ncbi/CF-random_software)  

* **Modeling Active-State Conformations of G-Protein-Coupled Receptors Using AlphaFold2 via Template Bias and Explicit Protein Constrains** [2025]  
Luca Chiesa, Dina Khasanova, and Esther Kellenberger.  
[J. Chem. Inf. Model. (2025)](https://doi.org/10.1021/acs.jcim.5c00489) | [code](https://github.com/LIT-CCM-lab/LIT-AlphaFold)  

* **AlphaFold prediction of structural ensembles of disordered proteins** [2025]  
Brotzakis, Z.F., Zhang, S., Murtada, M.H. et al.  
[Nat Commun 16, 1632 (2025)](https://doi.org/10.1038/s41467-025-56572-9) | [code](https://github.com/vendruscolo-lab/AlphaFold-IDP)  

* **Hidden Structural States of Proteins Revealed by Conformer Selection with AlphaFold-NMR** [2025]  
Yuanpeng J. Huang, Theresa A. Ramelot, Laura E. Spaman, Naohiro Kobayashi, Gaetano T. Montelione.  
[bioRxiv (2025)](https://doi.org/10.1101/2024.06.26.600902) | [code](https://github.rpi.edu/RPIBioinformatics/AlphaFold-NMR)  

* **AFflecto: A web server to generate conformational ensembles of flexible proteins from AlphaFold models** [2025]  
Pajkos, Mátyás, Ilinka Clerc, Christophe Zanon, Pau Bernadó, and Juan Cortés.  
[Journal of Molecular Biology (2025)](https://doi.org/10.1016/j.jmb.2025.169003) | [web](https://moma.laas.fr/applications/AFflecto)  

* **Characterizing the Conformational States of G Protein Coupled Receptors Generated with AlphaFold** [2025]  
Garima Chib, Parisa Mollaei, Amir Barati Farimani.  
[arXiv:2502.17628(2025)](hhttps://arxiv.org/abs/2502.17628) | [code](https://github.com/garimachib01/GPCR_AlphaFold)  

* **Gradations in protein dynamics captured by experimental NMR are not well represented by AlphaFold2 models and other computational metrics** [2025]  
Gavalda-Garcia, Jose, Bhawna Dixit, Adrián Díaz, An Ghysels, and Wim Vranken.   
[Journal of Molecular Biology 437.2 (2025)](https://doi.org/10.1016/j.jmb.2024.168900) | [code](https://zenodo.org/records/14192555)  

* **Modeling Protein Conformations by Guiding AlphaFold2 with Distance Distributions. Application to Double Electron Electron Resonance (DEER) Spectroscopy** [2024]  
Tianqi Wu, Richard A. Stein, Te-Yu Kao, Benjamin Brown, Hassane S. Mchaourab   
[bioRxiv. (2024)](https://doi.org/10.1101/2024.10.30.621127)  

* **AlphaFold-Multimer accurately captures interactions and dynamics of intrinsically disordered protein regions** [2024]  
Alireza Omidi, Mads Harder Møller, Nawar Malhis, and Jörg Gsponer.   
[bioRxiv. (2024)](https://doi.org/10.1073/pnas.2406407121) | [code](https://github.com/alirezaomidi/AFM-IDR)  

* **Harnessing AlphaFold to reveal hERG channel conformational state secrets** [2024]  
Khoa Ngo, Pei-Chi Yang, Vladimir Yarov-Yarovoy, Colleen E. Clancy, Igor Vorobyov.   
[bioRxiv. (2024)](https://doi.org/10.1101/2024.01.27.577468)  

* **AlphaFold2's training set powers its predictions of fold-switched conformations** [2024]  
Joseph W. Schafer, Lauren Porter.   
[bioRxiv. (2024)](https://doi.org/10.1101/2024.10.11.617857) | [data](https://github.com/porterll/CFold_AF2)  

* **AlphaFold2 Predicts Alternative Conformation Populations in Green Fluorescent Protein Variants** [2024]  
Núñez-Franco, Reyes, M. Milagros Muriel-Olaya, Gonzalo Jiménez-Osés, and Francesca Peccati.   
[J. Chem. Inf. Model. (2024)](https://doi.org/10.1021/acs.jcim.4c01388) | [data](https://zenodo.org/records/13133811)  

* **AlphaFold Ensemble Competition Screens Enable Peptide Binder Design with Single-Residue Sensitivity** [2024]  
Vosbein, Pernille, Paula Paredes Vergara, Danny T. Huang, and Andrew R. Thomson.   
[ACS Chemical Biology (2024)](https://doi.org/10.1016/j.str.2024.09.001)  

* **Assessing AF2’s ability to predict structural ensembles of proteins** [2024]  
Riccabona, Jakob R., Fabian C. Spoendlin, Anna-Lena M. Fischer, Johannes R. Loeffler, Patrick K. Quoika, Timothy P. Jenkins, James A. Ferguson et al.   
[Structure (2024)](https://doi.org/10.1016/j.str.2024.09.001)  

* **AlphaFold with conformational sampling reveals the structural landscape of homorepeats** [2024]  
Bonet, David Fernandez et al.   
[Structure (2024)](https://doi.org/10.1016/j.str.2024.08.016) | [code](	https://doi.org/10.5281/zenodo.13255318)  

* **Structure prediction of alternative protein conformations** [2024]  
Bryant, P., Noé, F.   
[Nat Commun 15, 7328 (2024)](https://doi.org/10.1038/s41467-024-51507-2) | [code](https://github.com/patrickbryant1/Cfold)  

* **AlphaFold predictions of fold-switched conformations are driven by structure memorization** [2024]  
Chakravarty, D., Schafer, J.W., Chen, E.A. et al.   
[Nat Commun 15, 7296 (2024)](https://doi.org/10.1038/s41467-024-51801-z) | [code](https://github.com/ncbi/AF2_benchmark)  

* **Predicting protein conformational motions using energetic frustration analysis and AlphaFold2** [2024]  
Xingyue Guan  and Qian-Yuan Tang  and Weitong Ren  and Mingchen Chen  and Wei Wang  and Peter G. Wolynes  and Wenfei Li.   
[Proceedings of the National Academy of Sciences (2024)](https://doi.org/10.1073/pnas.2410662121)  

* **Leveraging Machine Learning and AlphaFold2 Steering to Discover State-Specific Inhibitors Across the Kinome** [2024]  
Francesco Trozzi, Oanh Tran, Carmen Al Masri, Shu-Hang Lin, Balaguru Ravikumar, Rayees Rahman.   
[bioRxiv (2024)](https://doi.org/10.1101/2024.08.16.608358)  

* **A resource for comparing AF-Cluster and other AlphaFold2 sampling methods** [2024]  
Hannah K Wayment-Steele, Sergey Ovchinnikov, Lucy Colwell, Dorothee Kern.   
[bioRxiv (2024)](https://doi.org/10.1101/2024.07.29.605333)  

* **Integration of AlphaFold with Molecular Dynamics for Efficient Conformational Sampling of Transporter Protein NarK** [2024]  
Ohnuki, Jun, and Kei-ichi Okazaki.   
[The Journal of Physical Chemistry B (2024)](https://doi.org/10.1021/acs.jpcb.4c02726)  

* **AFsample2: Predicting multiple conformations and ensembles with AlphaFold2** [2024]  
Yogesh Kalakoti, Björn Wallner.   
[bioRxiv (2024)](https://doi.org/10.1101/2024.05.28.596195) | [code](https://github.com/iamysk/AFsample2/)  

* **Prediction of Conformational Ensembles and Structural Effects of State-Switching Allosteric Mutants in the Protein Kinases Using Comparative Analysis of AlphaFold2 Adaptations with Sequence Masking and Shallow Subsampling** [2024]  
Nishank Raisinghani, Mohammed Alshahrani, Grace Gupta, Hao Tian, Sian Xiao, Peng Tao, Gennady Verkhivker.   
[bioRxiv (2024)](https://doi.org/10.1101/2024.05.17.594786) | [code](https://zenodo.org/records/11204773)

* **Empowering AlphaFold2 for protein conformation selective drug discovery with AlphaFold2-RAVE** [2024]  
Xinyu Gu, Akashnathan Aranganathan, Pratyush Tiwary.   
[	arXiv:2404.07102 (2024)](https://arxiv.org/abs/2404.07102)  

* **High-throughput prediction of protein conformational distributions with subsampled AlphaFold2** [2024]  
Monteiro da Silva, G., Cui, J.Y., Dalgarno, D.C. et al.   
[Nat Commun 15, 2464 (2024)](https://doi.org/10.1038/s41467-024-46715-9) | [code](https://github.com/GMdSilva/gms_natcomms_1705932980_data)

* **AlphaFold Meets Flow Matching for Generating Protein Ensembles** [2024]  
Jing, Bowen, Bonnie Berger, and Tommi Jaakkola.   
[arXiv:2402.04845 (2024)](https://arxiv.org/abs/2402.04845) | [code](https://github.com/bjing2016/alphaflow)

* **Predicting multiple conformations via sequence clustering and AlphaFold2** [2024]  
Wayment-Steele, H.K., Ojoawo, A., Otten, R. et al.   
[Nature 625, 832–839 (2024)](https://doi.org/10.1038/s41586-023-06832-9) | [code](https://github.com/HWaymentSteele/AF_Cluster)

* **AlphaFold2-RAVE: From Sequence to Boltzmann Ranking** [2023]  
Bodhi P. Vani, Akashnathan Aranganathan, Dedi Wang, and Pratyush Tiwary.   
[J. Chem. Theory Comput. (2023)](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00290)) | [code](https://github.com/tiwarylab/alphafold2rave)

* **Investigating the conformational landscape of AlphaFold2-predicted protein kinase structures** [2023]  
Carmen Al-Masri, Francesco Trozzi, Shu-Hang Lin, Oanh Tran, Navriti Sahni, Marcel Patek, Anna Cichonska, Balaguru Ravikumar, Rayees Rahman.   
[Bioinformatics Advances. (2023)](https://doi.org/10.1093/bioadv/vbad129)) | [code](https://github.com/Harmonic-Discovery/AF2-kinase-conformational-landscape)

* **Exploring the Druggable Conformational Space of Protein Kinases Using AI-Generated Structures** [2023]  
Herrington, Noah B., David Stein, Yan Chak Li, Gaurav Pandey, and Avner Schlessinger.   
[bioRxiv (2023)](https://doi.org/10.1101/2023.08.31.555779) | [code](https://github.com/schlessinger-lab/af2_kinase_conformations/)

* **Sampling alternative conformational states of transporters and receptors with AlphaFold2** [2022]  
Del Alamo, Diego, Davide Sala, Hassane S. Mchaourab, and Jens Meiler.   
[Elife 11 (2022)](https://elifesciences.org/articles/75751) | [code](https://github.com/delalamo/af2_conformations)









### Autoregressive-Based







* **Simultaneous Modeling of Protein Conformation and Dynamics via Autoregression** [2025]   
Yuning Shen, Lihao Wang, Huizhuo Yuan, Yan Wang, Bangji Yang, Quanquan Gu.  
[arXiv:2505.17478 (2025)](https://doi.org/10.48550/arXiv.2505.17478)  

* **Force-Free Molecular Dynamics Through Autoregressive Equivariant Networks** [2025]   
Fabian L. Thiemann, Thiago Reschützegger, Massimiliano Esposito, Tseden Taddese, Juan D. Olarte-Plata, Fausto Martelli.  
[arXiv:2503.23794 (2025)](https://www.arxiv.org/abs/2503.23794) |  [code](https://github.com/IBM/trajcast)  



### LSTM-based







* **Learning molecular dynamics with simple language model built upon long short-term memory neural network** [2020]  
Tsai, ST., Kuo, EJ. & Tiwary, P.   
[Nat Commun 11, 5115 (2020)](https://doi.org/10.1038/s41467-020-18959-8) | [code](https://github.com/tiwarylab/LSTM-predict-MD)







### Transformer-based



* **Accurate Prediction of the Kinetic Sequence of Physicochemical States Using Generative Artificial Intelligence** [2025]  
Bera, Palash and Mondal, Jagannath.   
[Chem. Sci. (2025)](https://doi.org/10.1039/D5SC00108K) | [code](https://github.com/palash892/gpt_state_generation)  

* **Exploring the conformational ensembles of protein-protein complex with transformer-based generative model** [2024]  
Wang, Jianmin, Xun Wang, Yanyi Chu, Chunyan Li, Xue Li, Xiangyu Meng, Yitian Fang, Kyoung Tai No, Jiashun Mao, and Xiangxiang Zeng.   
[J. Chem. Theory Comput. (2024)](https://doi.org/10.1021/acs.jctc.4c00255) | [bioRxiv (2024)](https://doi.org/10.1101/2024.02.24.581708) | [code](https://github.com/AspirinCode/AlphaPPImd)  

* **Data-Efficient Generation of Protein Conformational Ensembles with Backbone-to-Side-Chain Transformers** [2024]  
Chennakesavalu, Shriram, and Grant M. Rotskoff.   
[The Journal of Physical Chemistry B (2024)](https://pubs.acs.org/doi/full/10.1021/acs.jpcb.3c08195) | [code](https://github.com/rotskoff-group/transformer-backmapping)

* **Molecular dynamics without molecules: searching the conformational space of proteins with generative neural networks** [2022]  
Schwing, Gregory, Luigi L. Palese, Ariel Fernández, Loren Schwiebert, and Domenico L. Gatti.   
[arXiv:2206.04683 (2022)](https://arxiv.org/abs/2206.04683) | [code](https://github.com/dgattiwsu/MD_without_molecules)








### VAE-based

* **Reinforced molecular dynamics: Physics-infused generative machine learning model explores CRBN activation process** [2025]  
Talant Ruzmetov, Ta I Hung, Saisri Padmaja Jonnalagedda, Si-han Chen, Parisa Fasihianifard, Zhefeng Guo, Bir Bhanu, Chia-en A. Chang.  
[bioRxiv. (2025)](https://doi.org/10.1101/2025.02.12.638002)  

* **Sampling Conformational Ensembles of Highly Dynamic Proteins via Generative Deep Learning** [2024]  
Talant Ruzmetov, Ta I Hung, Saisri Padmaja Jonnalagedda, Si-han Chen, Parisa Fasihianifard, Zhefeng Guo, Bir Bhanu, Chia-en A. Chang.  
[bioRxiv. (2024)](https://doi.org/10.1101/2024.05.05.592587) | [data](https://github.com/chang-group/ICoN)  

* **Exploring Protein Conformational Changes Using a Large-Scale Biophysical Sampling Augmented Deep Learning Strategy** [2024]  
Hu, Yao, Hao Yang, Mingwei Li, Zhicheng Zhong, Yongqi Zhou, Fang Bai, and Qian Wang.   
[Advanced Science (2024)](https://doi.org/10.1002/advs.202400884) | [code](https://github.com/qwang897/PATHpre)  

* **Deciphering the Coevolutionary Dynamics of L2 β-Lactamases via Deep Learning** [2024]  
  Zhu, Yu, Jing Gu, Zhuoran Zhao, AW Edith Chan, Maria F. Mojica, Andrea M. Hujer, Robert A. Bonomo, and Shozeb Haider.  
  [J. Chem. Inf. Model. (2024)](https://doi.org/10.1021/acs.jcim.4c00189) | [data](https://zenodo.org/records/10500539)

* **Protein Ensemble Generation Through Variational Autoencoder Latent Space Sampling** [2024]  
Sanaa Mansoor, Minkyung Baek, Hahnbeom Park, Gyu Rie Lee, and David Baker.   
[J. Chem. Theory Comput. (2024)](https://pubs.acs.org/doi/10.1021/acs.jctc.3c01057)  

* **Emerging Frontiers in Conformational Exploration of Disordered Proteins: Integrating Autoencoder and Molecular Simulations** [2024]  
Sanaa Mansoor, Minkyung Baek, Hahnbeom Park, Gyu Rie Lee, and David Baker.  
[ACS Chem. Neurosci. (2024)](https://doi.org/10.1021/acschemneuro.4c00670)  

* **Phanto-IDP: compact model for precise intrinsically disordered protein backbone generation and enhanced sampling** [2024]  
  Junjie Zhu, Zhengxin Li, Haowei Tong, Zhouyu Lu, Ningjie Zhang, Ting Wei and Hai-Feng Chen.  
  [Briefings in Bioinformatics. (2024)](https://academic.oup.com/bib/article/25/1/bbad429/7453435) | [code](https://github.com/Junjie-Zhu/Phanto-IDP)

* **Enhancing Conformational Sampling for Intrinsically Disordered and Ordered Proteins by Variational Auotencoder** [2023]  
  JunJie Zhu, NingJie Zhang, Ting Wei and Hai-Feng Chen.  
  [International Journal of Molecular Sciences. (2023)](https://www.mdpi.com/1422-0067/24/8/6896) | [code](https://github.com/Junjie-Zhu/VAE)

* **Encoding the Space of Protein-protein Binding Interfaces by Artificial Intelligence** [2023]  
  Su, Zhaoqian, Kalyani Dhusia, and Yinghao Wu.  
  [bioRxiv (2023)](https://doi.org/10.1101/2023.09.08.556812)  

* **Artificial intelligence guided conformational mining of intrinsically disordered proteins** [2022]  
Gupta, A., Dey, S., Hicks, A. et al.   
[Commun Biol 5, 610 (2022)](https://doi.org/10.1038/s42003-022-03562-y) | [code](https://github.com/aaayushg/generative_IDPs)

* **LAST: Latent Space-Assisted Adaptive Sampling for Protein Trajectories** [2022]  
Tian, Hao, Xi Jiang, Sian Xiao, Hunter La Force, Eric C. Larson, and Peng Tao   
[J. Chem. Inf. Model. (2022)](https://doi.org/10.1021/acs.jcim.2c01213) | [code](https://github.com/smu-tao-group/LAST)

* **Molecular dynamics without molecules: searching the conformational space of proteins with generative neural networks** [2022]  
Schwing, Gregory, Luigi L. Palese, Ariel Fernández, Loren Schwiebert, and Domenico L. Gatti.   
[arXiv:2206.04683 (2022)](https://arxiv.org/abs/2206.04683) | [code](https://github.com/dgattiwsu/MD_without_molecules)

* **ProGAE: A Geometric Autoencoder-based Generative Model for Disentangling Protein Conformational Space** [2021]  
Tatro, Norman Joseph, Payel Das, Pin-Yu Chen, Vijil Chenthamarakshan, and Rongjie Lai.   
[ICLR (2022)](https://openreview.net/forum?id=LxhlyKH6VP)

* **Explore protein conformational space with variational autoencoder** [2021]  
  Tian, Hao, Xi Jiang, Francesco Trozzi, Sian Xiao, Eric C. Larson, and Peng Tao.   
  [Frontiers in molecular biosciences 8 (2021)](https://www.frontiersin.org/articles/10.3389/fmolb.2021.781635/full) | [code](https://github.com/smu-tao-group/protein-VAE)








### GAN-based






* **DeepPath: Overcoming data scarcity for protein transition pathway prediction using physics-based deep learning** [2025]  
Yui Tik Pang, Katie M. Kuo, Lixinhao Yang, James C. Gumbart.   
[bioRxiv. (2025)](https://doi.org/10.1101/2025.02.27.640693)   

* **Direct generation of protein conformational ensembles via machine learning** [2023]  
Janson, G., Valdes-Garcia, G., Heo, L. et al.   
[Nat Commun 14, 774 (2023)](https://doi.org/10.1038/s41467-023-36443-x) | [code](https://github.com/feiglab/idpgan)  

* **Molecular dynamics without molecules: searching the conformational space of proteins with generative neural networks** [2022]  
Schwing, Gregory, Luigi L. Palese, Ariel Fernández, Loren Schwiebert, and Domenico L. Gatti.   
[arXiv:2206.04683 (2022)](https://arxiv.org/abs/2206.04683) | [code](https://github.com/dgattiwsu/MD_without_molecules)  









### Flow-based



* **Operator Forces For Coarse-Grained Molecular Dynamics** [2025]  
Leon Klein, Atharva Kelkar, Aleksander Durumeric, Yaoyi Chen, Frank Noé.   
[arXiv:2506.19628 (2025)](https://doi.org/10.48550/arXiv.2506.19628) | [data](https://github.com/noegroup/OperatorForces4CG)  

* **Enhanced Sampling, Public Dataset and Generative Model for Drug-Protein Dissociation Dynamics** [2025]  
Maodong Li, Jiying Zhang, Bin Feng, Wenqi Zeng, Dechin Chen, Zhijun Pan, Yu Li, Zijing Liu, Yi Isaac Yang.   
[arXiv:2504.18367 (2025)](https://arxiv.org/abs/2504.18367) | [data](https://huggingface.co/SZBL-IDEA)  

* **P2DFlow: A Protein Ensemble Generative Model with SE(3) Flow Matching** [2024]  
Yaowei Jin, Qi Huang, Ziyang Song, Mingyue Zheng, Dan Teng, Qian Shi.   
[	arXiv:2411.17196 (2024)](https://arxiv.org/abs/2411.17196) | [code](https://github.com/BLEACH366/P2DFlow)  

* **Generative Modeling of Molecular Dynamics Trajectories** [2024]  
Jing, Bowen, Hannes Stark, Tommi Jaakkola, and Bonnie Berger.   
[ICML'24 Workshop ML for Life and Material Science: From Theory to Industry Applications (2024)](https://openreview.net/forum?id=LbwM4VCDUU) | [code](https://github.com/bjing2016/mdgen)  

* **Frame-to-Frame Coarse-grained Molecular Dynamics with SE (3) Guided Flow Matching** [2024]  
Li, Shaoning, Yusong Wang, Mingyu Li, Jian Zhang, Bin Shao, Nanning Zheng, and Jian Tang   
[arXiv:2405.00751 (2024)](https://arxiv.org/abs/2405.00751)  

* **AlphaFold Meets Flow Matching for Generating Protein Ensembles** [2024]  
Jing, Bowen, Bonnie Berger, and Tommi Jaakkola.   
[arXiv:2402.04845 (2024)](https://arxiv.org/abs/2402.04845) | [code](https://github.com/bjing2016/alphaflow)  






### Diffusion-based





* **Generative Modeling of Full-Atom Protein Conformations using Latent Diffusion on Graph Embeddings** [2025]  
Aditya Sengar, Ali Hariri, Daniel Probst, Patrick Barth, Pierre Vandergheynst.  
[arXiv:2506.17064 (2025)](https://www.arxiv.org/abs/2506.17064) | [code](https://github.com/adityasengar/LD-FPG)  

* **Simultaneous Modeling of Protein Conformation and Dynamics via Autoregression** [2025]   
Yuning Shen, Lihao Wang, Huizhuo Yuan, Yan Wang, Bangji Yang, Quanquan Gu.  
[arXiv:2505.17478 (2025)](https://doi.org/10.48550/arXiv.2505.17478)  

* **Protein Conformation Generation via Force-Guided SE(3) Diffusion Models** [2024]  
Yan Wang, Lihao Wang, Yuning Shen, Yiqun Wang, Huizhuo Yuan, Yue Wu, Quanquan Gu.   
[ICML 2024 (2024)](https://arxiv.org/abs/2403.14088) | [code](https://github.com/bytedance/ConfDiff)  

* **Generative modeling of protein ensembles guided by crystallographic electron densities** [2024]  
Sai Advaith Maddipatla, Nadav Bojan Sellam, Sanketh Vedula, Ailie Marx, Alex Bronstein.   
[arXiv:2412.13223(2024)](https://arxiv.org/abs/2412.13223)  

* **Deep learning of protein energy landscape and conformational dynamics from experimental structures in PDB** [2024]  
Yike Tang, Mendi Yu, Ganggang Bai, Xinjun Li, Yanyan Xu, Buyong Ma.   
[bioRxiv (2024)](https://doi.org/10.1101/2024.06.27.600251)  

* **4D Diffusion for Dynamic Protein Structure Prediction with Reference Guided Motion Alignment** [2024]  
Cheng, Kaihui, Ce Liu, Qingkun Su, Jun Wang, Liwei Zhang, Yining Tang, Yao Yao, Siyu Zhu, and Yuan Qi.   
[arXiv:2408.12419 (2024)](https://arxiv.org/abs/2408.12419)  

* **Generating Multi-state Conformations of P-type ATPases with a Diffusion Model** [2024]  
Jingtian Xu, Yong Wang.   
[bioRxiv (2024)](https://doi.org/10.1101/2024.08.07.607107) | [code](https://github.com/yongwangCPH/papers/tree/main/2024/PtypeATPaseGeneration)  

* **Transferable deep generative modeling of intrinsically disordered protein conformations** [2024]  
Abdin, O., Kim, P.M.   
[PLOS Computational Biology 20.5 (2024)](https://doi.org/10.1371/journal.pcbi.1012144) | [code](https://github.com/giacomo-janson/idpsam)  

* **Direct conformational sampling from peptide energy landscapes through hypernetwork-conditioned diffusion** [2024]  
Janson, Giacomo, and Michael Feig.   
[Nat Mach Intell 6, 775–786 (2024)](https://doi.org/10.1038/s42256-024-00860-4) | [code](https://gitlab.com/oabdin/pepflow)

* **Accurate Conformation Sampling via Protein Structural Diffusion** [2024]  
Fan, Jiahao, Ziyao Li, Eric Alcaide, Guolin Ke, Huaqing Huang, and Weinan E.   
[J. Chem. Inf. Model. (2024)](https://doi.org/10.1021/acs.jcim.4c00928) | [bioRxiv (2024)](https://doi.org/10.1101/2024.05.20.594916) | [code](https://github.com/PKUfjh/Uni-Fold/tree/ufconf)    

* **Accurate and Efficient Structural Ensemble Generation of Macrocyclic Peptides using Internal Coordinate Diffusion** [2023]  
Grambow, Colin A., Hayley Weir, Nathaniel Diamant, Alex Tseng, Tommaso Biancalani, Gabriele Scalia and Kangway V Chuang.   
[	arXiv:2305.19800 (2023)](https://arxiv.org/abs/2305.19800) | [code](https://github.com/genentech/ringer)  






### Score-based






* **Str2str: A score-based framework for zero-shot protein conformation sampling** [2024]  
Lu, Jiarui, Bozitao Zhong, Zuobai Zhang, and Jian Tang.   
[ICLR (2024)](https://openreview.net/forum?id=C4BikKsgmK) | [code](https://github.com/lujiarui/Str2Str)

* **Score-based enhanced sampling for protein molecular dynamics** [2023]  
Lu, Jiarui, Bozitao Zhong, and Jian Tang.   
[arXiv:2306.03117 (2023)](https://arxiv.org/abs/2306.03117) | [code](https://github.com/lujiarui/Str2Str)






### Energy-based




* **Aligning Protein Conformation Ensemble Generation with Physical Feedback** [2025]  
Du, Yilun, Joshua Meier, Jerry Ma, Rob Fergus, and Alexander Rives.   
[ICML 2025 (2025)](https://openreview.net/forum?id=Asr955jcuZ) | [arXiv:2505.24203 (2025)](https://doi.org/10.48550/arXiv.2505.24203) | [code](https://github.com/lujiarui/eba)  

* **Energy-based models for atomic-resolution protein conformations** [2020]  
Du, Yilun, Joshua Meier, Jerry Ma, Rob Fergus, and Alexander Rives.   
[ICLR (2020)](https://openreview.net/forum?id=S1e_9xrFvS) | [code](https://github.com/facebookresearch/protein-ebm)  











### Bayesian-based


* **BaNDyT: Bayesian Network Modeling of Molecular Dynamics Trajectories** [2025]  
Mukhaleva, Elizaveta, Babgen Manookian, Hanyu Chen, Indira R. Sivaraj, Ning Ma, Wenyuan Wei, Konstancja Urbaniak et al.   
[J. Chem. Inf. Model. (2025)](https://doi.org/10.1021/acs.jcim.4c01981) | [code](https://github.com/bandyt-group/bandyt)  

* **Enabling Population Protein Dynamics Through Bayesian Modeling** [2024]  
Sylvain Lehmann, Jérôme Vialaret, Audrey Gabelle, Luc Bauchet, Jean-Philippe Villemin, Christophe Hirtz, Jacques Colinge.   
[Bioinformatics (2024)](https://doi.org/10.1093/bioinformatics/btae484)  

* **Deep Boosted Molecular Dynamics (DBMD): Accelerating molecular simulations with Gaussian boost potentials generated using probabilistic Bayesian deep neural network** [2023]  
Do, Hung N., and Yinglong Miao.   
[bioRxiv(2023)](https://doi.org/10.1101/2023.03.25.534210) | [code](https://github.com/MiaoLab20/DBMD/)  

* **Deep Generative Models of Protein Structure Uncover Distant Relationships Across a Continuous Fold Space** [2023]  
Draizen, Eli J., Stella Veretnik, Cameron Mura, and Philip E. Bourne.   
[bioRxiv(2023)](https://doi.org/10.1101/2022.07.29.501943) | [code](https://github.com/bouralab/DeepUrfold)  






### Active Learning-based



* **Active Learning of the Conformational Ensemble of Proteins Using Maximum Entropy VAMPNets** [2023]  
Kleiman, Diego E., and Diwakar Shukla.   
[J. Chem. Theory Comput. (2023)](https://doi.org/10.1021/acs.jctc.3c00040) | [code](https://github.com/ShuklaGroup/MaxEntVAMPNet)  







### GNN-based




* **A graph neural network-state predictive information bottleneck (GNN-SPIB) approach for learning molecular thermodynamics and kinetics** [2024]   
 Zou, Ziyue, Dedi Wang, and Pratyush Tiwary.  
  [Digital Discovery (2024)](https://doi.org/10.1039/D4DD00315B) | [code](https://github.com/connorzzou/PLUMED-NEST)  

* **Graph theory approaches for molecular dynamics simulations** [2024]   
 Patel AC, Sinha S, Palermo G.  
  [Quarterly Reviews of Biophysics. (2024)](https://doi.org/10.1017/S0033583524000143)  

* **EquiJump: Protein Dynamics Simulation via SO(3)-Equivariant Stochastic Interpolants** [2024]  
Allan dos Santos Costa, Ilan Mitnikov, Franco Pellegrini, Ameya Daigavane, Mario Geiger, Zhonglin Cao, Karsten Kreis, Tess Smidt, Emine Kucukbenli, Joseph Jacobson.   
[	arXiv:2410.09667 (2024)](https://arxiv.org/abs/2410.09667)  

* **AbFlex: Predicting the conformational flexibility of antibody CDRs** [2024]  
Spoendlin, Fabian C., Wing Ki Wong, Guy Georges, Alexander Bujotzek, and Charlotte Deane.   
[ICML'24 Workshop ML for Life and Material Science: From Theory to Industry Applications (2024)](https://openreview.net/forum?id=or4tArwd5a) | [code](https://openreview.net/forum?id=or4tArwd5a)

* **RevGraphVAMP: A protein molecular simulation analysis model combining graph convolutional neural networks and physical constraints** [2024]  
Huang, Ying, Huiling Zhang, Zhenli Lin, Yanjie Wei, and Wenhui Xi.   
[bioRxiv (2024)](https://doi.org/10.1101/2024.03.11.584426) | [code](https://github.com/DS00HY/RevGraphVamp)  




### LLM-MD


* **MD-LLM-1: A Large Language Model for Molecular Dynamics** [2025]  
Mhd Hussein Murtada, Z. Faidon Brotzakis, Michele Vendruscolo.   
[arXiv:2508.03709 (2025)](https://doi.org/10.48550/arXiv.2508.03709)  

* **Structure Language Models for Protein Conformation Generation** [2024]  
Jiarui Lu, Xiaoyin Chen, Stephen Zhewen Lu, Chence Shi, Hongyu Guo, Yoshua Bengio, Jian Tang.   
[arXiv:2410.18403 (2024)](https://arxiv.org/abs/2410.18403) | [code](https://github.com/lujiarui/esmdiff)  

* **SeaMoon: Prediction of molecular motions based on language models** [2024]  
Valentin Lombard, Dan Timsit, Sergei Grudinin, Elodie Laine.   
[bioRxiv. (2024)](https://doi.org/10.1101/2024.09.23.614585) | [code](https://github.com/PhyloSofS-Team//seamoon)  

* **Molecular simulation with an LLM-agent** [2024]  
MD-Agent is a LLM-agent based toolset for Molecular Dynamics.  
[code](https://github.com/ur-whitelab/md-agent)









## Molecular conformational ensembles by methods




### Small molecule conformational ensembles





* **Template-Guided 3D Molecular Pose Generation via Flow Matching and Differentiable Optimization** [2025]  
Bergues, Noémie, Arthur Carré, Paul Join-Lambert, Brice Hoffmann, Arnaud Blondel, and Hamza Tajmouati..  
[arXiv:2506.06305 (2025)](https://doi.org/10.48550/arXiv.2506.06305) |  [code](https://zenodo.org/records/15395813)  

* **DihedralsDiff: A Diffusion Conformation Generation Model That Unifies Local and Global Molecular Structures** [2025]  
Jianhui Xiao, Zheng Zheng, and Hao Liu.  
[J. Chem. Inf. Model. (2025)](https://doi.org/10.1021/acs.jcim.5c00367) |  [code](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/JNGTDF)  

* **WGFormer: An SE(3)-Transformer Driven by Wasserstein Gradient Flows for Molecular Ground-State Conformation Prediction** [2025]  
Fanmeng Wang, Minjie Cheng, Hongteng Xu.  
[ICML (2025)](https://doi.org/10.48550/arXiv.2410.09795) |  [code](https://github.com/FanmengWang/WGFormer)  

* **The pucke.rs toolkit to facilitate sampling the conformational space of biomolecular monomers** [2025]  
Rihon, J., Reynders, S., Bernardes Pinheiro, V. et al.  
[J Cheminform 17, 53 (2025)](https://doi.org/10.1186/s13321-025-00977-7) |  [code](https://github.com/jrihon/puckepy)  

* **Challenges and opportunities for machine learning potentials in transition path sampling: alanine dipeptide and azobenzene studies** [2025]  
Fedik, Nikita and Li, Wei and Lubbers, Nicholas and Nebgen, Benjamin and Tretiak, Sergei and Li, Ying Wai.  
[Digital Discovery (2025)](https://doi.org/10.1039/D4DD00265B) |  [code](https://github.com/nikitafedik/ml_tps_si)  

* **Diffusion-based generative AI for exploring transition states from 2D molecular graphs** [2024]   
Kim, S., Woo, J. & Kim, W.Y.   
[Nat Commun 15, 341 (2024)](https://doi.org/10.1038/s41467-023-44629-6) |  [code](https://github.com/seonghann/tsdiff)  

* **Physics-informed generative model for drug-like molecule conformers** [2024]   
David C. Williams, Neil Imana.   
[arXiv:2403.07925. (2024)](https://arxiv.org/abs/2403.07925v1) |  [code](https://github.com/nobiastx/diffusion-conformer)   

* **ET-Flow: Equivariant Flow-Matching for Molecular Conformer Generation** [2024]   
Majdi Hassan, Nikhil Shenoy, Jungyoon Lee, Hannes Stark, Stephan Thaler, Dominique Beaini.   
[NeurIPS 2024 (2024)](https://arxiv.org/abs/2410.22388) |  [code](https://github.com/shenoynikhil/ETFlow)  

* **COSMIC: Molecular Conformation Space Modeling in Internal Coordinates with an Adversarial Framework** [2024]   
Kuznetsov, Maksim, Fedor Ryabov, Roman Schutski, Rim Shayakhmetov, Yen-Chu Lin, Alex Aliper, and Daniil Polykovskiy.   
[J. Chem. Inf. Model. (2024)](https://doi.org/10.1021/acs.jcim.3c00989) |  [code](https://github.com/insilicomedicine/COSMIC)  

* **Leveraging 2D Molecular Graph Pretraining for Improved 3D Conformer Generation with Graph Neural Networks** [2024]   
Jiang, Runxuan, Tarun Gogineni, Joshua Kammeraad, Yifei He, Ambuj Tewari, and Paul M. Zimmerman.   
[Computers & Chemical Engineering (2024)](https://doi.org/10.1016/j.compchemeng.2024.108622) |  [code](https://github.com/m1k2zoo/2D-3DConformerGNN) 

* **DynamicsDiffusion: Generating and Rare Event Sampling of Molecular Dynamic Trajectories Using Diffusion Models** [2023]   
Petersen, Magnus, Gemma Roig, and Roberto Covino.   
[NeurIPS 2023 AI4Science  (2023)](https://openreview.net/forum?id=pwYCCq4xAf) 

* **Generating Molecular Conformer Fields** [2023]   
Yuyang Wang, Ahmed Elhag, Navdeep Jaitly, Joshua Susskind, Miguel Bautista.   
[NeurIPS 2023 Generative AI and Biology (GenBio) Workshop (2023)]https://openreview.net/forum?id=Od1KtMeAYo) 

* **On Accelerating Diffusion-based Molecular Conformation Generation in SE(3)-invariant Space** [2023]   
Zhou, Z., Liu, R. and Yu, T.   
[arXiv:2310.04915 (2023))](https://arxiv.org/abs/2310.04915) 

* **Molecular Conformation Generation via Shifting Scores** [2023]   
Zhou, Zihan, Ruiying Liu, Chaolong Ying, Ruimao Zhang, and Tianshu Yu.   
[arXiv:2309.09985 (2023)](https://arxiv.org/abs/2309.09985) 

* **EC-Conf: An Ultra-fast Diffusion Model for Molecular Conformation Generation with Equivariant Consistency** [2023]   
Fan, Zhiguang, Yuedong Yang, Mingyuan Xu, and Hongming Chen.   
[arXiv:2308.00237 (2023)](https://arxiv.org/abs/2308.00237) 

* **Prediction of Molecular Conformation Using Deep Generative Neural Networks** [2023]   
Xu, Congsheng, Yi Lu, Xiaomei Deng, and Peiyuan Yu.   
[Chinese Journal of Chemistry(2023)](https://doi.org/10.1002/cjoc.202300269)

* **Learning Over Molecular Conformer Ensembles: Datasets and Benchmarks** [2023]   
Zhu, Yanqiao, Jeehyun Hwang, Keir Adams, Zhen Liu, Bozhao Nan, Brock Stenfors, Yuanqi Du et al.   
[NeurIPS 2023 AI for Science Workshop. 2023 (2023)](https://openreview.net/forum?id=kFiMXnLH9x) |  [code](https://github.com/SXKDZ/MARCEL) 

* **Deep-Learning-Assisted Enhanced Sampling for Exploring Molecular Conformational Changes** [2023]   
Haohao Fu, Han Liu, Jingya Xing, Tong Zhao, Xueguang Shao, and Wensheng Cai.   
[J. Phys. Chem. B (2023)](https://doi.org/10.1021/acs.jpcb.3c05284) 

* **Torsional diffusion for molecular conformer generation** [2022]   
Jing, Bowen, Gabriele Corso, Jeffrey Chang, Regina Barzilay, and Tommi Jaakkola.   
[NeurIPS. (2022)](https://proceedings.neurips.cc/paper_files/paper/2022/hash/994545b2308bbbbc97e3e687ea9e464f-Abstract-Conference.html) |  [code](https://github.com/gcorso/torsional-diffusionf) 

* **GeoDiff: A Geometric Diffusion Model for Molecular Conformation Generation** [2022]   
Xu, Minkai, Lantao Yu, Yang Song, Chence Shi, Stefano Ermon, and Jian Tang.   
[International Conference on Learning Representations. (2022)](https://openreview.net/forum?id=PzcvxEMzvQC) |  [code](https://github.com/MinkaiXu/GeoDiff) 

* **Conformer-RL: A deep reinforcement learning library for conformer generation** [2022]   
Jiang, Runxuan, Tarun Gogineni, Joshua Kammeraad, Yifei He, Ambuj Tewari, and Paul M. Zimmerman.   
[Journal of Computational Chemistry 43.27 (2022)](https://doi.org/10.1002/jcc.26984) |  [code](https://github.com/ZimmermanGroup/conformer-rl) 

* **Energy-inspired molecular conformation optimization** [2022]   
Guan, Jiaqi, Wesley Wei Qian, Wei-Ying Ma, Jianzhu Ma, and Jian Peng.   
[International Conference on Learning Representations. (2022)](https://openreview.net/forum?id=7QfLW-XZTl) |  [code](https://github.com/guanjq/confopt_officialf) 

* **An End-to-End Framework for Molecular Conformation Generation via Bilevel Programming** [2021]   
Xu, Minkai, Wujie Wang, Shitong Luo, Chence Shi, Yoshua Bengio, Rafael Gomez-Bombarelli, and Jian Tang.   
[International Conference on Machine Learning. PMLR (2021)](http://proceedings.mlr.press/v139/xu21f.html) |  [code](https://github.com/MinkaiXu/ConfVAE-ICML21)  








### RNA conformational ensembles


* **DynaRNA: Dynamic RNA Conformation Ensemble Generation with Diffusion Model** [2025]  
Zhengxin Li, Junjie Zhu, Xiaokun Hong, Zhuoqi Zheng, Taeyoung Cui, Yutong Sun, Ting Wei, Haifeng Chen.   
[bioRxiv. (2025)](https://doi.org/10.1101/2025.05.22.655453)  

* **Determining structures of RNA conformers using AFM and deep neural networks** [2025]  
Boccalini, Matteo, Yelyzaveta Berezovska, Giovanni Bussi, Matteo Paloni, and Alessandro Barducci.  
[Proceedings of the National Academy of Sciences 122.15 (2025)](https://doi.org/10.1073/pnas.2425261122) | [code](https://doi.org/10.5281/zenodo.14203857)  

* **Determining structures of RNA conformers using AFM and deep neural networks** [2024]  
Degenhardt, M.F.S., Degenhardt, H.F., Bhandari, Y.R. et al.   
[Nature (2024)](https://doi.org/10.1038/s41586-024-07559-x) | [code](https://zenodo.org/records/10637777)  

* **On the Power and Challenges of Atomistic Molecular Dynamics to Investigate RNA Molecules** [2024]  
Muscat, Stefano, Gianfranco Martino, Jacopo Manigrasso, Marco Marcia, and Marco De Vivo.   
[J. Chem. Theory Comput. (2024)](https://doi.org/10.1021/acs.jctc.4c00773)  

* **Conformational ensembles of RNA oligonucleotides from integrating NMR and molecular simulations** [2018]  
Bottaro, S., Bussi, G., Kennedy, S.D., Turner, D.H. and Lindorff-Larsen, K.   
[Science advances 4.5 (2018)](https://doi.org/10.1038/s41597-024-03698-y) | [code](https://github.com/sbottaro/rr) | [data](https://github.com/sbottaro/tetranucleotides_data)  







### Peptide conformational ensembles





* **Memory kernel minimization-based neural networks for discovering slow collective variables of biomolecular dynamics** [2025]  
Liu, B., Cao, S., Boysen, J.G. et al.  
[Nat Comput Sci (2025)](https://doi.org/10.1038/s43588-025-00815-8) | [code](https://github.com/markovmodel/mdshare)  

* **Scoring Conformational Metastability of Macrocyclic Peptides with Binding Pose Metadynamics** [2025]  
Ryan Dykstra and Dan Sindhikara.  
[J. Chem. Inf. Model. (2025)](https://doi.org/10.1021/acs.jcim.4c01408)  

* **CREMP: Conformer-rotamer ensembles of macrocyclic peptides for machine learning** [2024]  
Grambow, C.A., Weir, H., Cunningham, C.N. et al.   
[Sci Data 11, 859 (2024)](https://doi.org/10.1038/s41597-024-03698-y) | [code](https://github.com/Genentech/cremp)  

* **Direct conformational sampling from peptide energy landscapes through hypernetwork-conditioned diffusion** [2024]  
Abdin, O., Kim, P.M.   
[Nat Mach Intell 6, 775–786 (2024)](https://doi.org/10.1038/s42256-024-00860-4) | [code](https://gitlab.com/oabdin/pepflow)  

* **Accurate and Efficient Structural Ensemble Generation of Macrocyclic Peptides using Internal Coordinate Diffusion** [2023]  
Grambow, Colin A., Hayley Weir, Nathaniel Diamant, Alex Tseng, Tommaso Biancalani, Gabriele Scalia and Kangway V Chuang.   
[	arXiv:2305.19800 (2023)](https://arxiv.org/abs/2305.19800) | [code](https://github.com/genentech/ringer)  









### Protein conformational ensembles








* **Applied causality to infer protein dynamics and kinetics** [2025]  
Akashnathan Aranganathan, Eric R. Beyerle.  
[arXiv:2508.12060 (2025)](https://doi.org/10.48550/arXiv.2508.12060)  

* **af2rave: Protein Ensemble Generation with Physics-Based Sampling** [2025]  
Teng D, Meraz VJ, Aranganathan A, Gu X, Tiwary P.  
[Digital Discovery (2025)](https://doi.org/10.1039/D5DD00201J) | [ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-q3mwr) | [code](https://github.com/tiwarylab/af2rave)  

* **Generative Modeling of Full-Atom Protein Conformations using Latent Diffusion on Graph Embeddings** [2025]  
Aditya Sengar, Ali Hariri, Daniel Probst, Patrick Barth, Pierre Vandergheynst.  
[arXiv:2506.17064 (2025)](https://www.arxiv.org/abs/2506.17064) | [code](https://github.com/adityasengar/LD-FPG)  

* **Memory kernel minimization-based neural networks for discovering slow collective variables of biomolecular dynamics** [2025]  
Liu, B., Cao, S., Boysen, J.G. et al.  
[Nat Comput Sci (2025)](https://doi.org/10.1038/s43588-025-00815-8) | [code](https://github.com/markovmodel/mdshare)  

* **Aligning Protein Conformation Ensemble Generation with Physical Feedback** [2025]  
Du, Yilun, Joshua Meier, Jerry Ma, Rob Fergus, and Alexander Rives.   
[ICML 2025 (2025)](https://openreview.net/forum?id=Asr955jcuZ) | [arXiv:2505.24203 (2025)](https://doi.org/10.48550/arXiv.2505.24203) | [code](https://github.com/lujiarui/eba)  

* **Simultaneous Modeling of Protein Conformation and Dynamics via Autoregression** [2025]   
Yuning Shen, Lihao Wang, Huizhuo Yuan, Yan Wang, Bangji Yang, Quanquan Gu.  
[arXiv:2505.17478 (2025)](https://doi.org/10.48550/arXiv.2505.17478)  

* **Caver Web 2.0: analysis of tunnels and ligand transport in dynamic ensembles of proteins** [2025]  
Sérgio M Marques, Simeon Borko, Ondrej Vavra, Jan Dvorsky, Petr Kohout, Petr Kabourek, Lukas Hejtmanek, Jiri Damborsky, David Bednar.   
[Nucleic Acids Research (2025)](https://doi.org/10.1093/nar/gkaf399) | [web](https://loschmidt.chemi.muni.cz/caverweb)  

* **GōMartini 3: From large conformational changes in proteins to environmental bias corrections** [2025]  
Souza, P.C.T., Borges-Araújo, L., Brasnett, C. et al.   
[Nat Commun 16, 4051 (2025)](https://doi.org/10.1038/s41467-025-58719-0) | [code](https://github.com/marrink-lab/vermouth-martinize)  

* **Learning Biophysical Dynamics with Protein Language Models** [2025]  
Chao Hou, Haiqing Zhao, Yufeng Shen.  
[bioRxiv. (2025)](https://doi.org/10.1101/2024.10.11.617911)  

* **Emerging Frontiers in Conformational Exploration of Disordered Proteins: Integrating Autoencoder and Molecular Simulations** [2024]  
Sanaa Mansoor, Minkyung Baek, Hahnbeom Park, Gyu Rie Lee, and David Baker.  
[ACS Chem. Neurosci. (2024)](https://doi.org/10.1021/acschemneuro.4c00670)  

* **Accurate Prediction of the Kinetic Sequence of Physicochemical States Using Generative Artificial Intelligence** [2025]  
Bera, Palash and Mondal, Jagannath.   
[Chem. Sci. (2025)](https://doi.org/10.1039/D5SC00108K) | [code](https://github.com/palash892/gpt_state_generation)  

* **Towards a Unified Framework for Determining Conformational Ensembles of Disordered Proteins** [2025]  
Hamidreza Ghafouri and Pavel Kadeřávek and Ana M Melo. et al.  
[arXiv:2504.03590 (2025)](https://arxiv.org/abs/2504.03590)  

* **Protein Conformation Generation via Force-Guided SE(3) Diffusion Models** [2024]  
Yan Wang, Lihao Wang, Yuning Shen, Yiqun Wang, Huizhuo Yuan, Yue Wu, Quanquan Gu.   
[ICML 2024 (2024)](https://arxiv.org/abs/2403.14088) | [code](https://github.com/bytedance/ConfDiff)  

* **AlphaFold prediction of structural ensembles of disordered proteins** [2025]  
Brotzakis, Z.F., Zhang, S., Murtada, M.H. et al.  
[Nat Commun 16, 1632 (2025)](https://doi.org/10.1038/s41467-025-56572-9) | [code](https://github.com/vendruscolo-lab/AlphaFold-IDP)  

* **Unsupervised Learning of Progress Coordinates during Weighted Ensemble Simulations: Application to NTL9 Protein Folding** [2025]  
Jeremy M. G. Leung, Nicolas C. Frazee, Alexander Brace, Anthony T. Bogetti, Arvind Ramanathan, and Lillian T. Chong.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.4c01136) | [code](https://github.com/westpa/DL-enhancedWE)  

* **Generalizable Protein Dynamics in Serine-Threonine Kinases: Physics is the key** [2025]  
Soumendranath Bhakat, Shray Vats, Andreas Mardt, Alexei Degterev.  
[bioRxiv (2025)](https://doi.org/10.1101/2025.03.06.641878) | [code](https://github.com/svats73/mdml)  

* **Hidden Structural States of Proteins Revealed by Conformer Selection with AlphaFold-NMR** [2025]  
Yuanpeng J. Huang, Theresa A. Ramelot, Laura E. Spaman, Naohiro Kobayashi, Gaetano T. Montelione.  
[bioRxiv (2025)](https://doi.org/10.1101/2024.06.26.600902) | [code](https://github.rpi.edu/RPIBioinformatics/AlphaFold-NMR)  

* **Insights into phosphorylation-induced influences on conformations and inhibitor binding of CDK6 through GaMD trajectory-based deep learning** [2025]  
Zhao, Lu and Wang, Jian and Yang, Wanchun and Zhang, Canqing and Zhang, Weiwei and Chen, Jianzhong.   
[Phys. Chem. Chem. Phys. (2025)](http://dx.doi.org/10.1039/D4CP04579C) | [code](https://github.com/Peach92/CDK6-pho-GaMD-trajectory-based-deep-learning)  

* **Gradations in protein dynamics captured by experimental NMR are not well represented by AlphaFold2 models and other computational metrics** [2025]  
Gavalda-Garcia, Jose, Bhawna Dixit, Adrián Díaz, An Ghysels, and Wim Vranken.   
[Journal of Molecular Biology 437.2 (2025)](https://doi.org/10.1016/j.jmb.2024.168900) | [code](https://zenodo.org/records/14192555)  

* **Generative modeling of protein ensembles guided by crystallographic electron densities** [2024]  
Sai Advaith Maddipatla, Nadav Bojan Sellam, Sanketh Vedula, Ailie Marx, Alex Bronstein.   
[arXiv:2412.13223(2024)](https://arxiv.org/abs/2412.13223)  

* **Scalable emulation of protein equilibrium ensembles with generative deep learning** [2024]  
Sarah Lewis, Tim Hempel, José Jiménez Luna, Michael Gastegger, Yu Xie, Andrew Y. K. Foong, Victor García Satorras, Osama Abdin, Bastiaan S. Veeling, Iryna Zaporozhets, Yaoyi Chen, Soojung Yang, Arne Schneuing, Jigyasa Nigam, Federico Barbero, Vincent Stimper, Andrew Campbell, Jason Yim, Marten Lienen, Yu Shi, Shuxin Zheng, Hannes Schulz, Usman Munir, Cecilia Clementi, Frank Noé.   
[bioRxiv. (2024)](https://doi.org/10.1101/2024.12.05.626885)  

* **P2DFlow: A Protein Ensemble Generative Model with SE(3) Flow Matching** [2024]  
Yaowei Jin, Qi Huang, Ziyang Song, Mingyue Zheng, Dan Teng, Qian Shi.   
[arXiv:2411.17196 (2024)](https://arxiv.org/abs/2411.17196) | [code](https://github.com/BLEACH366/P2DFlow)  

* **Fast Sampling of Protein Conformational Dynamics** [2024]  
Michael A. Sauer, Souvik Mondal, Brandon Neff, Sthitadhi Maiti, Matthias Heyden.   
[	arXiv:2411.08154 (2024)](https://arxiv.org/abs/2411.08154)  

* **Sampling Conformational Ensembles of Highly Dynamic Proteins via Generative Deep Learning** [2024]  
Talant Ruzmetov, Ta I Hung, Saisri Padmaja Jonnalagedda, Si-han Chen, Parisa Fasihianifard, Zhefeng Guo, Bir Bhanu, Chia-en A. Chang.   
[bioRxiv. (2024)](https://doi.org/10.1101/2024.05.05.592587) | [data](https://github.com/chang-group/ICoN)  

* **AlphaFold2's training set powers its predictions of fold-switched conformations** [2024]  
Joseph W. Schafer, Lauren Porter.   
[bioRxiv. (2024)](https://doi.org/10.1101/2024.10.11.617857) | [data](https://github.com/porterll/CFold_AF2)  

* **Exploring Protein Conformational Changes Using a Large-Scale Biophysical Sampling Augmented Deep Learning Strategy** [2024]  
Hu, Yao, Hao Yang, Mingwei Li, Zhicheng Zhong, Yongqi Zhou, Fang Bai, and Qian Wang.   
[Advanced Science (2024)](https://doi.org/10.1002/advs.202400884) | [code](https://github.com/qwang897/PATHpre)  

* **AMARO: All Heavy-Atom Transferable Neural Network Potentials of Protein Thermodynamics** [2024]  
Mirarchi, Antonio, Raul P. Pelaez, Guillem Simeon, and Gianni De Fabritiis.  
[arXiv:2409.17852 (2024)](https://arxiv.org/abs/2409.17852) | [code](https://github.com/compsciencelab/amaro)  

* **Conformations of KRAS4B Affected by Its Partner Binding and G12C Mutation: Insights from GaMD Trajectory-Image Transformation-Based Deep Learning** [2024]  
Chen, Jianzhong, Jian Wang, Wanchun Yang, Lu Zhao, and Guodong Hu.   
[J. Chem. Inf. Model. (2024)](https://doi.org/10.1021/acs.jcim.4c01174) | [code](https://github.com/jianzhong70/Visualiztion-for-MD-trajectory-analysis)  

* **Assessing AF2’s ability to predict structural ensembles of proteins** [2024]  
Riccabona, Jakob R., Fabian C. Spoendlin, Anna-Lena M. Fischer, Johannes R. Loeffler, Patrick K. Quoika, Timothy P. Jenkins, James A. Ferguson et al.   
[Structure (2024)](https://doi.org/10.1016/j.str.2024.09.001)  

* **Identifying protein conformational states in the Protein Data Bank: Toward unlocking the potential of integrative dynamics studies** [2024]  
Ellaway, J. I., Anyango, S., Nair, S., Zaki, H. A., Nadzirin, N., Powell, H. R., ... & Velankar, S.   
[Structural Dynamics (2024)](https://doi.org/10.1063/4.0000251)   

* **AlphaFold with conformational sampling reveals the structural landscape of homorepeats** [2024]  
Bonet, David Fernandez et al.   
[Structure (2024)](https://doi.org/10.1016/j.str.2024.08.016) | [code](	https://doi.org/10.5281/zenodo.13255318)  

* **Structure prediction of alternative protein conformations** [2024]  
Bryant, P., Noé, F.   
[Nat Commun 15, 7328 (2024)](https://doi.org/10.1038/s41467-024-51507-2) | [code](https://github.com/patrickbryant1/Cfold)  

* **Deep learning guided design of dynamic proteins** [2024]  
Amy B. Guo, Deniz Akpinaroglu, Mark J.S. Kelly, Tanja Kortemme.   
[bioRxiv. (2024)](https://doi.org/10.1101/2024.07.17.603962)  

* **AlphaFold predictions of fold-switched conformations are driven by structure memorization** [2024]  
Chakravarty, D., Schafer, J.W., Chen, E.A. et al.   
[Nat Commun 15, 7296 (2024)](https://doi.org/10.1038/s41467-024-51801-z) | [code](https://github.com/ncbi/AF2_benchmark)  

* **4D Diffusion for Dynamic Protein Structure Prediction with Reference Guided Motion Alignment** [2024]  
Cheng, Kaihui, Ce Liu, Qingkun Su, Jun Wang, Liwei Zhang, Yining Tang, Yao Yao, Siyu Zhu, and Yuan Qi.   
[arXiv:2408.12419 (2024)](https://arxiv.org/abs/2408.12419)  

* **Predicting protein conformational motions using energetic frustration analysis and AlphaFold2** [2024]  
Xingyue Guan  and Qian-Yuan Tang  and Weitong Ren  and Mingchen Chen  and Wei Wang  and Peter G. Wolynes  and Wenfei Li.   
[Proceedings of the National Academy of Sciences (2024)](https://doi.org/10.1073/pnas.2410662121)  

* **A resource for comparing AF-Cluster and other AlphaFold2 sampling methods** [2024]  
Hannah K Wayment-Steele, Sergey Ovchinnikov, Lucy Colwell, Dorothee Kern.   
[bioRxiv (2024)](https://doi.org/10.1101/2024.07.29.605333) 

* **Integration of AlphaFold with Molecular Dynamics for Efficient Conformational Sampling of Transporter Protein NarK** [2024]  
Ohnuki, Jun, and Kei-ichi Okazaki.   
[The Journal of Physical Chemistry B (2024)](https://doi.org/10.1021/acs.jpcb.4c02726)  

* **Transferable deep generative modeling of intrinsically disordered protein conformations** [2024]  
Abdin, O., Kim, P.M.   
[PLOS Computational Biology 20.5 (2024)](https://doi.org/10.1371/journal.pcbi.1012144) | [code](https://github.com/giacomo-janson/idpsam)  

* **AFsample2: Predicting multiple conformations and ensembles with AlphaFold2** [2024]  
Yogesh Kalakoti, Björn Wallner.   
[bioRxiv (2024)](https://doi.org/10.1101/2024.05.28.596195) | [code](https://github.com/iamysk/AFsample2/)

* **Prediction of Conformational Ensembles and Structural Effects of State-Switching Allosteric Mutants in the Protein Kinases Using Comparative Analysis of AlphaFold2 Adaptations with Sequence Masking and Shallow Subsampling** [2024]  
Nishank Raisinghani, Mohammed Alshahrani, Grace Gupta, Hao Tian, Sian Xiao, Peng Tao, Gennady Verkhivker.   
[bioRxiv (2024)](https://doi.org/10.1101/2024.05.17.594786) | [code](https://zenodo.org/records/11204773)

* **Empowering AlphaFold2 for protein conformation selective drug discovery with AlphaFold2-RAVE** [2024]  
Xinyu Gu, Akashnathan Aranganathan, Pratyush Tiwary.   
[	arXiv:2404.07102 (2024)](https://arxiv.org/abs/2404.07102)  

* **High-throughput prediction of protein conformational distributions with subsampled AlphaFold2** [2024]  
Monteiro da Silva, G., Cui, J.Y., Dalgarno, D.C. et al.   
[Nat Commun 15, 2464 (2024)](https://doi.org/10.1038/s41467-024-46715-9) | [code](https://github.com/GMdSilva/gms_natcomms_1705932980_data)

* **AlphaFold Meets Flow Matching for Generating Protein Ensembles** [2024]  
Jing, Bowen, Bonnie Berger, and Tommi Jaakkola.   
[arXiv:2402.04845 (2024)](https://arxiv.org/abs/2402.04845) | [code](https://github.com/bjing2016/alphaflow)

* **Predicting multiple conformations via sequence clustering and AlphaFold2** [2024]  
Wayment-Steele, H.K., Ojoawo, A., Otten, R. et al.   
[Nature 625, 832–839 (2024)](https://doi.org/10.1038/s41586-023-06832-9) | [code](https://github.com/HWaymentSteele/AF_Cluster)

* **Data-Efficient Generation of Protein Conformational Ensembles with Backbone-to-Side-Chain Transformers** [2024]  
Chennakesavalu, Shriram, and Grant M. Rotskoff.   
[The Journal of Physical Chemistry B (2024)](https://pubs.acs.org/doi/full/10.1021/acs.jpcb.3c08195) | [code](https://github.com/rotskoff-group/transformer-backmapping)

* **Frame-to-Frame Coarse-grained Molecular Dynamics with SE (3) Guided Flow Matching** [2024]  
Li, Shaoning, Yusong Wang, Mingyu Li, Jian Zhang, Bin Shao, Nanning Zheng, and Jian Tang   
[arXiv:2405.00751 (2024)](https://arxiv.org/abs/2405.00751)  

* **AlphaFold Meets Flow Matching for Generating Protein Ensembles** [2024]  
Jing, Bowen, Bonnie Berger, and Tommi Jaakkola.   
[arXiv:2402.04845 (2024)](https://arxiv.org/abs/2402.04845) | [code](https://github.com/bjing2016/alphaflow)

* **Machine learning of force fields towards molecular dynamics simulations of proteins at DFT accuracy** [2024]  
Brunken, Christoph, Sebastien Boyer, Mustafa Omar, Bakary N'tji Diallo, Karim Beguir, Nicolas Lopez Carranza, and Oliver Bent.   
[ICLR 2024 Workshop on Generative and Experimental Perspectives for Biomolecular Design (2024)](https://openreview.net/forum?id=hrvvIOx7EM) | [code](https://github.com/ACEsuit/mace-jax)  

* **Deciphering the Coevolutionary Dynamics of L2 β-Lactamases via Deep Learning** [2024]  
  Zhu, Yu, Jing Gu, Zhuoran Zhao, AW Edith Chan, Maria F. Mojica, Andrea M. Hujer, Robert A. Bonomo, and Shozeb Haider.  
  [J. Chem. Inf. Model. (2024)](https://doi.org/10.1021/acs.jcim.4c00189) | [data](https://zenodo.org/records/10500539)

* **Protein Ensemble Generation Through Variational Autoencoder Latent Space Sampling** [2024]  
Sanaa Mansoor, Minkyung Baek, Hahnbeom Park, Gyu Rie Lee, and David Baker.   
[J. Chem. Theory Comput. (2024)](https://pubs.acs.org/doi/10.1021/acs.jctc.3c01057)

* **Phanto-IDP: compact model for precise intrinsically disordered protein backbone generation and enhanced sampling** [2024]  
  Junjie Zhu, Zhengxin Li, Haowei Tong, Zhouyu Lu, Ningjie Zhang, Ting Wei and Hai-Feng Chen.  
  [Briefings in Bioinformatics. (2024)](https://academic.oup.com/bib/article/25/1/bbad429/7453435) | [code](https://github.com/Junjie-Zhu/Phanto-IDP)

* **Str2str: A score-based framework for zero-shot protein conformation sampling** [2024]  
Lu, Jiarui, Bozitao Zhong, Zuobai Zhang, and Jian Tang.   
[ICLR (2024)](https://openreview.net/forum?id=C4BikKsgmK) | [code](https://github.com/lujiarui/Str2Str)

* **Enabling Population Protein Dynamics Through Bayesian Modeling** [2024]  
Sylvain Lehmann, Jérôme Vialaret, Audrey Gabelle, Luc Bauchet, Jean-Philippe Villemin, Christophe Hirtz, Jacques Colinge.   
[Bioinformatics (2024)](https://doi.org/10.1093/bioinformatics/btae484)

* **Deep Boosted Molecular Dynamics (DBMD): Accelerating molecular simulations with Gaussian boost potentials generated using probabilistic Bayesian deep neural network** [2023]  
Do, Hung N., and Yinglong Miao.   
[bioRxiv(2023)](https://doi.org/10.1101/2023.03.25.534210) | [code](https://github.com/MiaoLab20/DBMD/)

* **Deep Generative Models of Protein Structure Uncover Distant Relationships Across a Continuous Fold Space** [2023]  
Draizen, Eli J., Stella Veretnik, Cameron Mura, and Philip E. Bourne.   
[bioRxiv(2023)](https://doi.org/10.1101/2022.07.29.501943) | [code](https://github.com/bouralab/DeepUrfold)

* **Score-based enhanced sampling for protein molecular dynamics** [2023]  
Lu, Jiarui, Bozitao Zhong, and Jian Tang.   
[arXiv:2306.03117 (2023)](https://arxiv.org/abs/2306.03117) | [code](https://github.com/lujiarui/Str2Str)

* **AlphaFold2-RAVE: From Sequence to Boltzmann Ranking** [2023]  
Bodhi P. Vani, Akashnathan Aranganathan, Dedi Wang, and Pratyush Tiwary.   
[J. Chem. Theory Comput. (2023)](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00290)) | [code](https://github.com/tiwarylab/alphafold2rave)

* **Investigating the conformational landscape of AlphaFold2-predicted protein kinase structures** [2023]  
Carmen Al-Masri, Francesco Trozzi, Shu-Hang Lin, Oanh Tran, Navriti Sahni, Marcel Patek, Anna Cichonska, Balaguru Ravikumar, Rayees Rahman.   
[Bioinformatics Advances. (2023)](https://doi.org/10.1093/bioadv/vbad129)) | [code](https://github.com/Harmonic-Discovery/AF2-kinase-conformational-landscape)

* **Exploring the Druggable Conformational Space of Protein Kinases Using AI-Generated Structures** [2023]  
Herrington, Noah B., David Stein, Yan Chak Li, Gaurav Pandey, and Avner Schlessinger.   
[bioRxiv (2023)](https://doi.org/10.1101/2023.08.31.555779) | [code](https://github.com/schlessinger-lab/af2_kinase_conformations/)

* **Active Learning of the Conformational Ensemble of Proteins Using Maximum Entropy VAMPNets** [2023]  
Kleiman, Diego E., and Diwakar Shukla.   
[J. Chem. Theory Comput. (2023)](https://doi.org/10.1021/acs.jctc.3c00040) | [code](https://github.com/ShuklaGroup/MaxEntVAMPNet)

* **Direct generation of protein conformational ensembles via machine learning** [2023]  
Janson, G., Valdes-Garcia, G., Heo, L. et al.   
[Nat Commun 14, 774 (2023)](https://doi.org/10.1038/s41467-023-36443-x) | [code](https://github.com/feiglab/idpgan)

* **Enhancing Conformational Sampling for Intrinsically Disordered and Ordered Proteins by Variational Auotencoder** [2023]  
  JunJie Zhu, NingJie Zhang, Ting Wei and Hai-Feng Chen.  
  [International Journal of Molecular Sciences. (2023)](https://www.mdpi.com/1422-0067/24/8/6896) | [code](https://github.com/Junjie-Zhu/VAE)

* **Molecular dynamics without molecules: searching the conformational space of proteins with generative neural networks** [2022]  
Schwing, Gregory, Luigi L. Palese, Ariel Fernández, Loren Schwiebert, and Domenico L. Gatti.   
[arXiv:2206.04683 (2022)](https://arxiv.org/abs/2206.04683) | [code](https://github.com/dgattiwsu/MD_without_molecules)

* **Sampling alternative conformational states of transporters and receptors with AlphaFold2** [2022]  
Del Alamo, Diego, Davide Sala, Hassane S. Mchaourab, and Jens Meiler.   
[Elife 11 (2022)](https://elifesciences.org/articles/75751) | [code](https://github.com/delalamo/af2_conformations)  

* **Artificial intelligence guided conformational mining of intrinsically disordered proteins** [2022]  
Gupta, A., Dey, S., Hicks, A. et al.   
[Commun Biol 5, 610 (2022)](https://doi.org/10.1038/s42003-022-03562-y) | [code](https://github.com/aaayushg/generative_IDPs)

* **LAST: Latent Space-Assisted Adaptive Sampling for Protein Trajectories** [2022]  
Tian, Hao, Xi Jiang, Sian Xiao, Hunter La Force, Eric C. Larson, and Peng Tao   
[J. Chem. Inf. Model. (2022)](https://doi.org/10.1021/acs.jcim.2c01213) | [code](https://github.com/smu-tao-group/LAST)

* **Molecular dynamics without molecules: searching the conformational space of proteins with generative neural networks** [2022]  
Schwing, Gregory, Luigi L. Palese, Ariel Fernández, Loren Schwiebert, and Domenico L. Gatti.   
[arXiv:2206.04683 (2022)](https://arxiv.org/abs/2206.04683) | [code](https://github.com/dgattiwsu/MD_without_molecules)

* **ProGAE: A Geometric Autoencoder-based Generative Model for Disentangling Protein Conformational Space** [2021]  
Tatro, Norman Joseph, Payel Das, Pin-Yu Chen, Vijil Chenthamarakshan, and Rongjie Lai.   
[ICLR (2022)](https://openreview.net/forum?id=LxhlyKH6VP)

* **Explore protein conformational space with variational autoencoder** [2021]  
  Tian, Hao, Xi Jiang, Francesco Trozzi, Sian Xiao, Eric C. Larson, and Peng Tao.   
  [Frontiers in molecular biosciences 8 (2021)](https://www.frontiersin.org/articles/10.3389/fmolb.2021.781635/full) | [code](https://github.com/smu-tao-group/protein-VAE)

* **Energy-based models for atomic-resolution protein conformations** [2020]  
Du, Yilun, Joshua Meier, Jerry Ma, Rob Fergus, and Alexander Rives.   
[ICLR (2020)](https://openreview.net/forum?id=S1e_9xrFvS) | [code](https://github.com/facebookresearch/protein-ebm)






### Enzymes conformational ensembles


* **Generating Multi-state Conformations of P-type ATPases with a Diffusion Model** [2024]  
Jingtian Xu, Yong Wang.   
[bioRxiv (2024)](https://doi.org/10.1101/2024.08.07.607107) | [code](https://github.com/yongwangCPH/papers/tree/main/2024/PtypeATPaseGeneration)

* **Deciphering the Coevolutionary Dynamics of L2 β-Lactamases via Deep Learning** [2024]  
  Zhu, Yu, Jing Gu, Zhuoran Zhao, AW Edith Chan, Maria F. Mojica, Andrea M. Hujer, Robert A. Bonomo, and Shozeb Haider.  
  [J. Chem. Inf. Model. (2024)](https://doi.org/10.1021/acs.jcim.4c00189) | [data](https://zenodo.org/records/10500539)








### Antibody conformational ensembles




* **AbFlex: Predicting the conformational flexibility of antibody CDRs** [2024]  
Spoendlin, Fabian C., Wing Ki Wong, Guy Georges, Alexander Bujotzek, and Charlotte Deane.   
[ICML'24 Workshop ML for Life and Material Science: From Theory to Industry Applications (2024)](https://openreview.net/forum?id=or4tArwd5a) | [code](https://openreview.net/forum?id=or4tArwd5a)








### Ligand-Protein conformational ensembles



* **Caver Web 2.0: analysis of tunnels and ligand transport in dynamic ensembles of proteins** [2025]  
Sérgio M Marques, Simeon Borko, Ondrej Vavra, Jan Dvorsky, Petr Kohout, Petr Kabourek, Lukas Hejtmanek, Jiri Damborsky, David Bednar.   
[Nucleic Acids Research (2025)](https://doi.org/10.1093/nar/gkaf399) | [web](https://loschmidt.chemi.muni.cz/caverweb)  

* **Enhanced Sampling, Public Dataset and Generative Model for Drug-Protein Dissociation Dynamics** [2025]  
Maodong Li, Jiying Zhang, Bin Feng, Wenqi Zeng, Dechin Chen, Zhijun Pan, Yu Li, Zijing Liu, Yi Isaac Yang.   
[arXiv:2504.18367 (2025)](https://arxiv.org/abs/2504.18367) | [data](https://huggingface.co/SZBL-IDEA)  

* **Molecular dynamics-powered hierarchical geometric deep learning framework for protein-ligand interaction** [2025]  
Liu, Mingquan and Jin, Shuting and Lai, Houtim and Wang, Longyue and Wang, Jianmin and Cheng, Zhixiang and Zeng, Xiangxiang.  
[IEEE Transactions on Computational Biology and Bioinformatics. (2025)](https://doi.org/10.1109/TCBBIO.2025.3558959) | [code](https://github.com/lmqfly/Dynamics-PLI)  

* **Towards automated physics-based absolute drug residence time predictions** [2025]  
Smith Z, Branduardi D, Lupyan D, D’Arrigo G, Tiwary P, Wang L, et al.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-wg75c)  

* **Comparative Analysis of Quantum-Mechanical and standard Single-Structure Protein-Ligand Scoring Functions with MD-Based Free Energy Calculations** [2025]  
SJalaie M, Fanfrlík J, Pecina A, Lepšík M, Řezáč J.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-38lf5) | [code](https://github.com/Honza-R/Wang_dataset_SQM)  

* **A Multi-Grained Symmetric Differential Equation Model for Learning Protein-Ligand Binding Dynamics** [2024]  
Shengchao Liu, Weitao Du, Hannan Xu, Yanjing Li, Zhuoxinran Li, Vignesh Bhethanabotla, Divin Yan, Christian Borgs, Anima Anandkumar, Hongyu Guo, Jennifer Chayes.   
[arXiv:2401.15122 (2024)](https://arxiv.org/abs/2401.15122) | [code](https://github.com/chao1224/NeuralMD)  

* **Modeling protein-small molecule conformational ensembles with ChemNet** [2024]  
Ivan Anishchenko, Yakov Kipnis, Indrek Kalvet, Guangfeng Zhou, Rohith Krishna, Samuel J. Pellock, Anna Lauko, Gyu Rie Lee, Linna An, Justas Dauparas, Frank DiMaio, David Baker.   
[bioRxiv (2024)](https://doi.org/10.1101/2024.09.25.614868) | [code](https://github.com/baker-laboratory/PLACER)  

* **MISATO: machine learning dataset of protein–ligand complexes for structure-based drug discovery** [2024]  
Siebenmorgen, T., Menezes, F., Benassou, S. et al.   
[Nat Comput Sci 4, 367–378 (2024)](https://doi.org/10.1038/s43588-024-00627-2) | [code](https://github.com/t7morgen/misato-dataset)  

* **Enhancing Protein–Ligand Binding Affinity Predictions Using Neural Network Potentials** [2024]  
Sabanés Zariquiey, F., Galvelis, R., Gallicchio, E., Chodera, J.D., Markland, T.E. and De Fabritiis, G.   
[J. Chem. Inf. Model. (2024)](https://doi.org/10.1021/acs.jcim.3c02031) | [code](https://github.com/compsciencelab/ATM_benchmark/tree/main/ATM_With_NNPs)  

* **Assessment of molecular dynamics time series descriptors in protein-ligand affinity prediction** [2024]  
Poziemski, Jakub, Artur Yurkevych, and Pawel Siedlecki.   
[chemrxiv-2024-dxv36 (2024)](https://doi.org/10.26434/chemrxiv-2024-dxv36) | [code](https://github.com/JPoziemski/md_for_affinity_prediction)  

* **Pre-Training of Equivariant Graph Matching Networks with Conformation Flexibility for Drug Binding** [2022]  
Wu, Fang, Shuting Jin, Yinghui Jiang, Xurui Jin, Bowen Tang, Zhangming Niu, Xiangrong Liu, Qiang Zhang, Xiangxiang Zeng, and Stan Z. Li.   
[Advanced Science 9.33 (2022)](https://doi.org/10.1002/advs.202203796) | [code](https://github.com/smiles724/ProtMD)  






### PPI conformational ensembles







* **Memory kernel minimization-based neural networks for discovering slow collective variables of biomolecular dynamics** [2025]  
Liu, B., Cao, S., Boysen, J.G. et al.  
[Nat Comput Sci (2025)](https://doi.org/10.1038/s43588-025-00815-8) | [code](https://github.com/markovmodel/mdshare)  

* **Generalizable Protein Dynamics in Serine-Threonine Kinases: Physics is the key** [2025]  
Soumendranath Bhakat, Shray Vats, Andreas Mardt, Alexei Degterev.  
[bioRxiv (2025)](https://doi.org/10.1101/2025.03.06.641878) | [code](https://github.com/svats73/mdml)  

* **Scalable emulation of protein equilibrium ensembles with generative deep learning** [2024]  
Sarah Lewis, Tim Hempel, José Jiménez Luna, Michael Gastegger, Yu Xie, Andrew Y. K. Foong, Victor García Satorras, Osama Abdin, Bastiaan S. Veeling, Iryna Zaporozhets, Yaoyi Chen, Soojung Yang, Arne Schneuing, Jigyasa Nigam, Federico Barbero, Vincent Stimper, Andrew Campbell, Jason Yim, Marten Lienen, Yu Shi, Shuxin Zheng, Hannes Schulz, Usman Munir, Cecilia Clementi, Frank Noé.   
[bioRxiv. (2024)](https://doi.org/10.1101/2024.12.05.626885)  

* **Computational screening of the effects of mutations on protein-protein off-rates and dissociation mechanisms by τRAMD** [2024]  
D’Arrigo, G., Kokh, D.B., Nunes-Alves, A. et al.   
[Commun Biol 7, 1159 (2024)](https://doi.org/10.1038/s42003-024-06880-5) | [code](https://github.com/HITS-MCM/tauRAMD)  

* **Quantifying conformational changes in the TCR:pMHC-I binding interface** [2024]  
Benjamin McMaster, Christopher Thorpe, Jamie Rossjohn, Charlotte M. Deane, Hashem Koohy.   
[bioRxiv (2024)](https://doi.org/10.1101/2024.08.13.607715) | [code](https://github.com/benjiemc/tcr-pmhc-interface-analysis)  

* **Exploring the conformational ensembles of protein-protein complex with transformer-based generative model** [2024]  
Wang, Jianmin, Xun Wang, Yanyi Chu, Chunyan Li, Xue Li, Xiangyu Meng, Yitian Fang, Kyoung Tai No, Jiashun Mao, and Xiangxiang Zeng.   
[J. Chem. Theory Comput. (2024)](https://doi.org/10.1021/acs.jctc.4c00255) | [bioRxiv (2024)](https://doi.org/10.1101/2024.02.24.581708) | [code](https://github.com/AspirinCode/AlphaPPImd)

* **Encoding the Space of Protein-protein Binding Interfaces by Artificial Intelligence** [2023]  
  Su, Zhaoqian, Kalyani Dhusia, and Yinghao Wu.  
  [bioRxiv (2023)](https://doi.org/10.1101/2023.09.08.556812)  

* **Unsupervised and supervised AI on molecular dynamics simulations reveals complex characteristics of HLA-A2-peptide immunogenicity** [2024]  
Jeffrey K Weber, Joseph A Morrone, Seung-gu Kang, Leili Zhang, Lijun Lang, Diego Chowell, Chirag Krishna, Tien Huynh, Prerana Parthasarathy, Binquan Luan, Tyler J Alban, Wendy D Cornell, Timothy A Chan.   
[Briefings in Bioinformatics (2024)](https://doi.org/10.1093/bib/bbad504) | [code](https://github.com/BiomedSciAI)  










### RNA-Peptide conformational ensembles






* **Enhanced Sampling Simulations of RNA-peptide Binding using Deep Learning Collective Variables** [2024]  
Nisha Kumari, Sonam Dhull, Tarak Karmakar.   
[bioRxiv (2024)](https://doi.org/10.1101/2024.08.01.606277)  







### Antibody-Protein conformational ensembles



* **Computational Mapping of Conformational Dynamics and Interaction Hotspots of Human VISTA with pH-Selective Antibodies** [2025]  
Norman Ly, Shubham Devesh Ramgoolam, and Aravindhan Ganesan.  
[Biochemistry (2025)](https://doi.org/10.1021/acs.biochem.4c00863)  

* **Using Short Molecular Dynamics Simulations to Determine the Important Features of Interactions in Antibody–Protein Complexes** [2024]  
A. Clay Richard, Robert J. Pantazes.   
[Proteins. (2024)](https://doi.org/10.1002/prot.26773)  








### Nucleic acid-Protein conformational ensembles




* **Communication pathway analysis within protein-nucleic acid complexes** [2025]  
Sneha Bheemireddy, Roy González-Alemán, Emmanuelle Bignon, Yasaman Karami.  
[bioRxiv. (2025)](https://doi.org/10.1101/2025.02.14.638259) | [code](https://github.com/yasamankarami/compass)  









### Nucleic acid-Ligand conformational ensembles





* **Machine learning-augmented molecular dynamics simulations (MD) reveal insights into the disconnect between affinity and activation of ZTP riboswitch ligands** [2025]  
Christopher Fullenkamp, Shams Mehdi, Christopher Jones, Logan Tenney, Patricio Pichling, Peri R. Prestwood, Adrian R. Ferré-D'Amaré, Pratyush Tiwary, John Schneekloth.  
[ Angew. Chem. Int. Ed. (2025)](https://doi.org/10.1002/anie.202505971)  






### Material ensembles



* **Accuracy of calculated elastic properties of inorganic materials: Density functional theory and machine-learned potentials** [2025]  
Milman, Victor, Alexander Perlov, Neil Spenley, and Björn Winkler.  
[Materialia (2025)](https://doi.org/10.1016/j.mtla.2025.102503)  

* **DREAMS: Density Functional Theory Based Research Engine for Agentic Materials Simulation** [2025]  
Ziqi Wang, Hongshuo Huang, Hancheng Zhao, Changwen Xu, Shang Zhu, Jan Janssen, Venkatasubramanian Viswanathan.  
[arXiv:2507.14267 (2025)](https://doi.org/10.48550/arXiv.2507.14267) | [code](https://github.com/BattModels/material_agent)  

* **Screening of Material Defects using Universal Machine-Learning Interatomic Potentials** [2025]  
PBerger, Ethan, Mohammad Bagheri, and Hannu-Pekka Komsa.  
[Small (2025)](https://doi.org/10.1002/smll.202503956) | [data](https://zenodo.org/records/15025795)  

* **chemtrain-deploy: A parallel and scalable framework for machine learning potentials in million-atom MD simulations** [2025]  
Paul Fuchs, Weilong Chen, Stephan Thaler, Julija Zavadlav.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.5c00996) | [arXiv:2506.04055 (2025)](https://doi.org/10.48550/arXiv.2506.04055) | [data](https://github.com/tummfm/chemtrain)  

* **Scalable Bayesian Optimization for High-Dimensional Coarse-Grained Model Parameterization** [2025]  
Carlos A. Martins Junior, Daniela A. Damasceno, Keat Yung Hue, Caetano R. Miranda, Erich A. Müller, Rodrigo A. Vargas-Hernández.  
[arXiv:2506.22533 (2025)](https://doi.org/10.48550/arXiv.2506.22533)  

* **PolyConf: Unlocking Polymer Conformation Generation through Hierarchical Generative Models** [2025]  
Fanmeng Wang, Wentao Guo, Qi Ou, Hongshuai Wang, Haitao Lin, Hongteng Xu, Zhifeng Gao.  
[ICML (2025)](https://arxiv.org/abs/2504.08859) | [code](https://github.com/FanmengWang/PolyConf)  

* **A predictive machine learning force-field framework for liquid electrolyte development** [2025]  
Gong, S., Zhang, Y., Mu, Z. et al.  
[Nat Mach Intell (2025)](https://doi.org/10.1038/s42256-025-01009-7) | [code](https://github.com/bytedance/bamboo)  

* **Machine learning-driven molecular dynamics unveils a bulk phase transformation driving ammonia synthesis on barium hydride** [2025]  
Tosello Gardini, A., Raucci, U. & Parrinello, M..  
[Nat Commun 16, 2475 (2025)](https://doi.org/10.1038/s41467-025-57688-8) | [code](https://doi.org/10.5281/zenodo.14789082)  

* **Advancing nonadiabatic molecular dynamics simulations in solids with E(3) equivariant deep neural hamiltonians** [2025]  
Zhang, C., Zhong, Y., Tao, ZG. et al.  
[Nat Commun 16, 2033 (2025)](https://doi.org/10.1038/s41467-025-57328-1) | [code](https://github.com/QuantumLab-ZY/HamGNN)  

* **End-To-End Learning of Classical Interatomic Potentials for Benchmarking Anion Polarization Effects in Lithium Polymer Electrolytes** [2025]  
Pablo A. Leon, Avni Singhal, Jurgis Ruza, Jeremiah A. Johnson, Yang Shao-Horn, and Rafael Gomez-Bombarelli.  
[Chem. Mater. (2025)](https://doi.org/10.1021/acs.chemmater.4c02529) | [code](https://github.com/learningmatter-mit/AutoBADDIE)   

* **General-purpose machine-learned potential for 16 elemental metals and their alloys** [2024]  
Song, K., Zhao, R., Liu, J. et al.  
[Nat Commun 15, 10208 (2024)](https://doi.org/10.1038/s41467-024-54554-x) | [code](https://github.com/brucefan1983/GPUMD)    

* **Quantum-accurate machine learning potentials for metal-organic frameworks using temperature driven active learning** [2024]  
Sharma, A., Sanvito, S.  
[npj Comput Mater 10, 237 (2024)](https://doi.org/10.1038/s41524-024-01427-y) | [code](https://github.com/asharma-ms/MOF_MLP_2024)  

* **Generative AI model trained by molecular dynamics for rapid mechanical design of architected graphene** [2024]  
Milad Masrouri, Kamalendu Paul, Zhao Qin.   
[Extreme Mechanics Letters (2024)](https://doi.org/10.1016/j.eml.2024.102230)  

* **Neural-network-based molecular dynamics simulations reveal that proton transport in water is doubly gated by sequential hydrogen-bond exchange** [2024]  
Gomez, A., Thompson, W.H. & Laage, D.   
[Nat. Chem. (2024)](https://doi.org/10.1038/s41557-024-01593-y) | [data](https://zenodo.org/records/11965260)  

* **Universal-neural-network-potential molecular dynamics for lithium metal and garnet-type solid electrolyte interface** [2024]  
Iwasaki, R., Tanibata, N., Takeda, H. et al.   
[Commun Mater 5, 148 (2024)](https://doi.org/10.1038/s43246-024-00595-0)  











