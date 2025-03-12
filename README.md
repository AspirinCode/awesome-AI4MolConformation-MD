[![License: GPL](https://img.shields.io/badge/License-GPL-yellow)](https://github.com/AspirinCode/awesome-AI4MolConformation-MD)

## awesome-AI4MolConformation-MD
List of **molecules ( small molecules, RNA, peptide, protein, enzymes, antibody, and PPIs) conformations** and **molecular dynamics (force fields)** using **generative artificial intelligence** and **deep learning**


![Protein Space and Conformations](https://github.com/AspirinCode/awesome-AI4MolConformation-MD/blob/main/figure/afpro.png)


**Updating ...**  


## Menu

  - [Deep Learning-molecular conformations](#deep-learning-molecular-conformations)

| Menu | Menu | Menu | Menu |
| ------ | :---------- | ------ | ------ |
| [Reviews](#reviews) | [Datasets and Package](#datasets-and-package) | [Molecular dynamics](#molecular-dynamics) | [Neural Molecular Force Fields](#neural-molecular-force-fields) |
| [MD Engines-Frameworks](#md-engines-frameworks) | [AI4MD Engines-Frameworks](#ai4md-engines-frameworks) | [MD Trajectory Processing-Analysis](#md-trajectory-processing-analysis) | [AI4MD](#ai4md) |
| [Neural Network Potentials](#neural-network-potentials) | [Free Energy Perturbation](#free-energy-perturbation) | [Ab Initio](#ab-initio) |  |
| [AlphaFold-based](#alphaFold-based) | [GNN-based](#gnn-based)  | [LSTM-based](#lstm-based) | [Transformer-based](#transformer-based) |
| [VAE-based](#vae-based) | [GAN-based](#gan-based) | [Flow-based](#flow-based) | [Diffusion-based](#diffusion-based) |
| [Score-Based](#score-Based) | [Energy-based](#energy-based) | [Bayesian-based](#bayesian-based) | [Active Learning-based](#active-learning-based) |
| [LLM-MD](#llm-md) |  |  |  |



  - [Molecular conformational ensembles by methods](#molecular-conformational-ensembles-by-methods)

| Menu | Menu | Menu |
| ------ | :---------- | ------ |
| [Small molecule conformational ensembles](#small-molecule-conformational-ensembles) | [RNA conformational ensembles](#rna-conformational-ensembles) | [Peptide conformational ensembles](#peptide-conformational-ensembles) |
| [Protein conformational ensembles](#protein-conformational-ensembles) | [Enzymes conformational ensembles](#enzymes-conformational-ensembles)  | [Antibody conformational ensembles](#antibody-conformational-ensembles) |
| [Ligand-Protein conformational ensembles](#ligand-protein-conformational-ensembles) | [RNA-Peptide conformational ensembles](#rna-peptide-conformational-ensembles) | [PPI conformational ensembles](#ppi-conformational-ensembles) |
| [Antibody-Protein conformational ensembles](#antibody-protein-conformational-ensembles) |  | [Material ensembles](#material-ensembles) |






## Reviews



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


* **OpenQDC: Open Quantum Data Commons** [2024]  
Cristian Gabellini, Nikhil Shenoy, Stephan Thaler, Semih Canturk, Daniel McNeela, Dominique Beaini, Michael Bronstein, Prudencio Tossou.   
[arXiv:2411.19629 (2024)](https://arxiv.org/abs/2411.19629) | [data](https://github.com/valence-labs/openQDC)  

* **Molecular Quantum Chemical Data Sets and Databases for Machine Learning Potentials** [2024]  
Antonio Mirarchi, Toni Giorgino, G. D. Fabritiis.   
[ChemRxiv. (2024)](https://doi.org/10.26434/chemrxiv-2024-w3ld0-v2) | [code](https://github.com/Arif-PhyChem/datasets_and_databases_4_MLPs)  

* **mdCATH: A Large-Scale MD Dataset for Data-Driven Computational Biophysics** [2024]  
Antonio Mirarchi, Toni Giorgino, G. D. Fabritiis.   
[	arXiv:2407.14794 (2024)](https://arxiv.org/abs/2407.14794) | [code](https://github.com/compsciencelab/mdCATH)  






### Package




**MMolearn**  
a Python package streamlining the design of generative models of biomolecular dynamics  

https://github.com/LumosBio/MolData   







## Molecular dynamics






### Neural Molecular Force Fields





* **Towards Fast, Specialized Machine Learning Force Fields: Distilling Foundation Models via Energy Hessians** [2025]  
Ishan Amin, Sanjeev Raja, Aditi Krishnapriyan.  
[arXiv:2501.09009 (2025)](https://arxiv.org/abs/2501.09009) | [code](https://github.com/ASK-Berkeley/MLFF-distill)  

* **BoostMD: Accelerated Molecular Sampling Leveraging ML Force Field Features** [2025]  
Zhaoxin Xie, Yanheng Li, Yijie Xia, Jun Zhang, Sihao Yuan, Cheng Fan, Yi Isaac Yang, and Yi Qin Gao.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.4c01449)  

* **Efficient and Precise Force Field Optimization for Biomolecules Using DPA-2** [2024]  
Junhan Chang, Duo Zhang, Yuqing Deng, Hongrui Lin, Zhirong Liu, Linfeng Zhang, Hang Zheng, Xinyan Wang.   
[arXiv:2406.09817 (2024)](https://arxiv.org/abs/2406.09817) | [code](https://github.com/deepmodeling/DeePDih)  

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

* **End-to-end differentiable construction of molecular mechanics force fields** [2022]  
Wang, Yuanqing, Josh Fass, Benjamin Kaminow, John E. Herr, Dominic Rufa, Ivy Zhang, Iván Pulido et al.   
[Chemical Science (2022)](https://doi.org/10.1039/D2SC02739A) | [code](https://github.com/openforcefield/openff-nagl)  


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

### MD Trajectory Processing-Analysis


* [MDAnalysis](https://www.mdanalysis.org/) - An object-oriented Python library to analyze trajectories from molecular dynamics (MD) simulations in many popular formats.  
* [MDTraj](http://mdtraj.org/) - A python library that allows users to manipulate molecular dynamics (MD) trajectories.  
* [PyTraj](https://amber-md.github.io/pytraj/) - A Python front-end package of the popular cpptraj program.  
* [CppTraj](https://github.com/Amber-MD/cpptraj) - Biomolecular simulation trajectory/data analysis.  
* [WEDAP](https://github.com/chonglab-pitt/wedap) - A Python Package for Streamlined Plotting of Molecular Simulation Data. 
* [Melodia](https://github.com/rwmontalvao/Melodia_py) - A Python library for protein structure analysis. 
* [MDANCE](https://github.com/mqcomplab/MDANCE) - A flexible n-ary clustering package that provides a set of tools for clustering Molecular Dynamics trajectories.  
* [PENSA](https://github.com/drorlab/pensa) - A collection of python methods for exploratory analysis and comparison of biomolecular conformational ensembles.  





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






## AI4MD



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

* **DeePMD-kit: A deep learning package for many-body potential energy representation and molecular dynamics** [2018]  
Wang, Han, Linfeng Zhang, Jiequn Han, and E. Weinan.   
[Computer Physics Communications 228 (2018)](https://doi.org/10.1016/j.cpc.2018.03.016) | [code](https://github.com/deepmodeling/deepmd-kit)  









### Neural Network Potentials





* **Free energy profiles for chemical reactions in solution from high-dimensional neural network potentials: The case of the Strecker synthesis** [2025]  
Alea Miako Tokita, Timothée Devergne, A. Marco Saitta, Jörg Behler.  
[arXiv:2503.05370 (2025)](https://arxiv.org/abs/2503.05370)  

* **Accurate Free Energy Calculation via Multiscale Simulations Driven by Hybrid Machine Learning and Molecular Mechanics Potentials** [2025]  
Wang X, Wu X, Brooks B, Wang J.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2024-zq975-v2)  

* **DeePMD-kit v3: A Multiple-Backend Framework for Machine Learning Potentials** [2025]  
Zeng, Jinzhe, Duo Zhang, Anyang Peng, Xiangyu Zhang, Sensen He, Yan Wang, Xinzijian Liu et al.  
[arXiv:2502.19161 (2025)](https://arxiv.org/abs/2502.19161) | [code](https://github.com/njzjz/benchmark-dpv3)  

* **ANI-1xBB: an ANI based reactive potential** [2025]  
Zhang S, Zubatyuk R, Yang Y, Roitberg A, Isayev O.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-m2nqq)  

* **Efficient Training of Neural Network Potentials for Chemical and Enzymatic Reactions by Continual Learning** [2025]  
Yao-Kun Lei, Kiyoshi Yagi, and Yuji Sugita.  
[J. Chem. Theory Comput. (2025)](https://doi.org/10.1021/acs.jctc.4c01393)   

* **Efficient equivariant model for machine learning interatomic potentials** [2025]  
Yang, Z., Wang, X., Li, Y. et al.  
[npj Comput Mater 11, 49 (2025)](https://doi.org/10.1038/s41524-025-01535-3) | [code](https://github.com/Shen-Group/E2GNN)   

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

* **DeepConf: Leveraging ANI-ML Potentials for Exploring Local Minima with A Focus on Bioactive Conformations** [2024]  
Tayfuroglu O, Zengin İN, Koca MS, Kocak A.  
[ChemRxiv. (2024)](https://doi.org/10.26434/chemrxiv-2024-xqbdp-v2) | [code](https://github.com/otayfuroglu/DeepConf)  

* **Universal Machine Learning Interatomic Potentials are Ready for Phonons** [2024]  
Antoine Loew, Dewen Sun, Hai-Chen Wang, Silvana Botti, Miguel A. L. Marques.  
[arXiv:2412.16551 (2024)](https://arxiv.org/abs/2412.16551) | [code](https://github.com/hyllios/utils)   

* **Neural Network Potentials for Enabling Advanced Small-Molecule Drug Discovery and Generative Design** [2024]  
Barnett, Simon, and John D. Chodera.  
[GEN Biotechnology 3.3 (2024)](https://doi.org/10.1089/genbio.2024.0011)    

* **General-purpose machine-learned potential for 16 elemental metals and their alloys** [2024]  
Song, K., Zhao, R., Liu, J. et al.  
[Nat Commun 15, 10208 (2024)](https://doi.org/10.1038/s41467-024-54554-x) | [code](https://github.com/brucefan1983/GPUMD)    

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

* **AIMNet2: A Neural Network Potential to Meet your Neutral, Charged, Organic, and Elemental-Organic Needs** [2024]  
Anstine, Dylan, Roman Zubatyuk, and Olexandr Isayev.   
[chemrxiv-2023-296ch-v2 (2024)](https://doi.org/10.26434/chemrxiv-2023-296ch-v2) | [code](https://github.com/isayevlab/aimnet2)  

* **NNP/MM: Accelerating Molecular Dynamics Simulations with Machine Learning Potentials and Molecular Mechanics** [2023]  
Galvelis, R., Varela-Rial, A., Doerr, S., Fino, R., Eastman, P., Markland, T.E., Chodera, J.D. and De Fabritiis, G.   
[J. Chem. Inf. Model. (2023)](https://doi.org/10.1021/acs.jcim.3c00773) | [code](https://github.com/openmm/nnpops)  

* **Towards universal neural network potential for material discovery applicable to arbitrary combination of 45 elements** [2022]  
Takamoto, S., Shinagawa, C., Motoki, D. et al.   
[Nat Commun 13, 2991 (2022)](https://doi.org/10.1038/s41467-022-30687-9) | [data](https://doi.org/10.6084/m9.figshare.19658538)  

* **Teaching a neural network to attach and detach electrons from molecules** [2021]  
Zubatyuk, R., Smith, J.S., Nebgen, B.T. et al.   
[Nat Commun 12, 4870 (2021)](https://doi.org/10.1038/s41467-021-24904-0)  | [code](https://github.com/isayevlab/aimnetnse)  

* **Four Generations of High-Dimensional Neural Network Potentials** [2021]  
Behler, Jorg.   
[Chemical Reviews 121.16 (2021)](https://doi.org/10.1021/acs.chemrev.0c00868)  

* **DP-GEN: A concurrent learning platform for the generation of reliable deep learning based potential energy models** [2020]  
Zhang, Yuzhi, Haidi Wang, Weijie Chen, Jinzhe Zeng, Linfeng Zhang, Han Wang, and E. Weinan.   
[Computer Physics Communications 253 (2020)](https://doi.org/10.1016/j.cpc.2020.107206)  | [code](https://github.com/deepmodeling/dpgen)  









### Free Energy Perturbation

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




* **Analytical ab initio hessian from a deep learning potential for transition state optimization** [2024]  
KYuan, E.CY., Kumar, A., Guan, X. et al.  
[Nat Commun 15, 8865 (2024)](https://doi.org/10.1038/s41467-024-52481-5) | [code](https://github.com/zadorlab/sella)  





## Deep Learning-molecular conformations


* **Molecular Simulations with a Pretrained Neural Network and Universal Pairwise Force Fields** [2024]  
Kabylda A, Frank JT, Dou SS, Khabibrakhmanov A, Sandonas LM, Unke OT, et al.  
[ChemRxiv. (2024)](https://doi.org/10.26434/chemrxiv-2024-bdfr0) | [code](https://github.com/general-molecular-simulations/so3lr)  

* **SpaiNN: Equivariant Message Passing for Excited-State Nonadiabatic Molecular Dynamics** [2024]  
Mausenberger, Sascha, Carolin Müller, Alexandre Tkatchenko, Philipp Marquetand, Leticia González, and Julia Westermayr.   
[Chemical Science (2024)](https://doi.org/10.1039/D4SC04164J) | [code](https://github.com/schnarc/SchNarc)  

* **GLOW: A Workflow Integrating Gaussian-Accelerated Molecular Dynamics and Deep Learning for Free Energy Profiling** [2022]  
Do, Hung N., Jinan Wang, Apurba Bhattarai, and Yinglong Miao.   
[J. Chem. Theory Comput. (2022)](https://doi.org/10.1021/acs.jctc.1c01055) | [code](https://github.com/MiaoLab20/GLOW)  



Pajkos, Mátyás, et al. "AFflecto: A web server to generate conformational ensembles of flexible proteins from AlphaFold models." Journal of Molecular Biology (2025): 169003.




### AlphaFold-based


* **Hidden Structural States of Proteins Revealed by Conformer Selection with AlphaFold-NMR** [2025]  
Yuanpeng J. Huang, Theresa A. Ramelot, Laura E. Spaman, Naohiro Kobayashi, Gaetano T. Montelione.  
[bioRxiv (2025)](https://doi.org/10.1101/2024.06.26.600902) | [code](https://github.rpi.edu/RPIBioinformatics/AlphaFold-NMR)  

* **AFflecto: A web server to generate conformational ensembles of flexible proteins from AlphaFold models** [2025]  
Pajkos, Mátyás, Ilinka Clerc, Christophe Zanon, Pau Bernadó, and Juan Cortés.  
[Journal of Molecular Biology (2025)](https://doi.org/10.1016/j.jmb.2025.169003) | [web](https://moma.laas.fr/applications/AFflecto)  

* **Characterizing the Conformational States of G Protein Coupled Receptors Generated with AlphaFold** [2025]  
Garima Chib, Parisa Mollaei, Amir Barati Farimani.  
[arXiv:2502.17628(2025)](hhttps://arxiv.org/abs/2502.17628) | [code](https://github.com/garimachib01/GPCR_AlphaFold)  

* **AlphaFold2-RAVE: Protein Ensemble Generation with Physics-Based Sampling** [2025]  
Teng D, Meraz VJ, Aranganathan A, Gu X, Tiwary P.  
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-q3mwr) | [code](https://github.com/tiwarylab/af2rave)  

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








### LSTM-based







* **Learning molecular dynamics with simple language model built upon long short-term memory neural network** [2020]  
Tsai, ST., Kuo, EJ. & Tiwary, P.   
[Nat Commun 11, 5115 (2020)](https://doi.org/10.1038/s41467-020-18959-8) | [code](https://github.com/tiwarylab/LSTM-predict-MD)







### Transformer-based





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



### LLM-MD




* **Structure Language Models for Protein Conformation Generation** [2024]  
Jiarui Lu, Xiaoyin Chen, Stephen Zhewen Lu, Chence Shi, Hongyu Guo, Yoshua Bengio, Jian Tang.   
[	arXiv:2410.18403 (2024)](https://arxiv.org/abs/2410.18403) | [code](https://github.com/lujiarui/esmdiff)  

* **SeaMoon: Prediction of molecular motions based on language models** [2024]  
Valentin Lombard, Dan Timsit, Sergei Grudinin, Elodie Laine.   
[bioRxiv. (2024)](https://doi.org/10.1101/2024.09.23.614585) | [code](https://github.com/PhyloSofS-Team//seamoon)  

* **Molecular simulation with an LLM-agent** [2024]  
MD-Agent is a LLM-agent based toolset for Molecular Dynamics.  
[code](https://github.com/ur-whitelab/md-agent)








## Molecular conformational ensembles by methods




### Small molecule conformational ensembles

* **Diffusion-based generative AI for exploring transition states from 2D molecular graphs** [2024]   
Kim, S., Woo, J. & Kim, W.Y.   
[Nat Commun 15, 341 (2024)](https://doi.org/10.1038/s41467-023-44629-6) |  [code](https://github.com/seonghann/tsdiff)

* **Physics-informed generative model for drug-like molecule conformers** [2024]   
David C. Williams, Neil Imana.   
[	arXiv:2403.07925. (2024)](https://arxiv.org/abs/2403.07925v1) |  [code](https://github.com/nobiastx/diffusion-conformer)   

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


* **Hidden Structural States of Proteins Revealed by Conformer Selection with AlphaFold-NMR** [2025]  
Yuanpeng J. Huang, Theresa A. Ramelot, Laura E. Spaman, Naohiro Kobayashi, Gaetano T. Montelione.  
[bioRxiv (2025)](https://doi.org/10.1101/2024.06.26.600902) | [code](https://github.rpi.edu/RPIBioinformatics/AlphaFold-NMR)  

* **Insights into phosphorylation-induced influences on conformations and inhibitor binding of CDK6 through GaMD trajectory-based deep learning** [2025]  
Zhao, Lu and Wang, Jian and Yang, Wanchun and Zhang, Canqing and Zhang, Weiwei and Chen, Jianzhong.   
[Phys. Chem. Chem. Phys. (2025)](http://dx.doi.org/10.1039/D4CP04579C) | [code](https://github.com/Peach92/CDK6-pho-GaMD-trajectory-based-deep-learning)  

* **AlphaFold2-RAVE: Protein Ensemble Generation with Physics-Based Sampling** [2025]  
Teng D, Meraz VJ, Aranganathan A, Gu X, Tiwary P.   
[ChemRxiv. (2025)](https://doi.org/10.26434/chemrxiv-2025-q3mwr) | [code](https://github.com/tiwarylab/af2rave)  

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





* **Using Short Molecular Dynamics Simulations to Determine the Important Features of Interactions in Antibody–Protein Complexes** [2024]  
A. Clay Richard, Robert J. Pantazes.   
[Proteins. (2024)](https://doi.org/10.1002/prot.26773)  




### Material ensembles



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











