## QSAR-hiv project 

Despite the significant progress in managing patients infected with HIV through the development of Highly Active Anti-Retroviral Therapy (HAART), major challenges and opportunities remain to be explored. Of particular interest, is the binding of glycoprotein 120 (gp120) to the primary cellular receptor Cluster of Differentiation 4 (CD4).  In this work we describe our  two phased computational process to identify useful compounds capable of binding to the gp120 protein for therapeutic purposes.  We identified 187 compounds from the literature that conform to active binding sites on these proteins and use these as training/test sets. The data in the form of quantitative structure-activity relationships (QSAR) is downloaded from the ZINC database and transformed using principal components analysis. In the first phase we developed a Radial Basis Function neural network model that identifies potential inhibitors from a virtual screen of a subset of the ZINC database. In the second phase we modelled the top performing compounds using the Discovery Studio docking and screening software. By employing this approach, we identified that those compounds with a LogP value of approx 2-4 performed well in the binding simulations while the lower scoring compounds do not bind well. 

This project is a collaboration between the University of Sunderland and the University of Newcastle.

If you find the data and R source code useful please cite our paper:

A. Hosny, M. Ashton, Y. Gong and K. McGarry, The development of a predictive model to identify potential HIV-1 attachment inhibitors, Accepted for publication in Computers in Biology and Medicine Journal, 2020.
