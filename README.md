This is a repository for my Diploma Thesis conducted in Aristotle University of Thessaloniki.
This thesis consists of two parts, namely URLLC Delay Minimization with FDMA and mMTC Sum Energy 
Consumption Minimization with NOMA.

Abstract

In recent years, the rapid increase of diverse wireless applications with various and competing service demands in terms of performance, bandwidth, reliability and latency has rendered the
existing «one-size-fits-all» 4G network design unable to meet the sheer need of supporting diverse quality of service (QoS) requirements. Furthermore, the traditional Cloud Computing architecture is unable to meet the performance requirements. In this context, radio access network
(RAN) slicing and Mobile Edge Computing (MEC) are two key enablers for 5G and beyond, particularly to empower low-latency services. This thesis investigates the coexistence of heterogeneous services on the same radio resource in an Edge Computing scenario. Services in 5G
are classified into three main categories, namely, Enhanced Mobile Broadband (eMBB), Massive Machine-Type Communications (mMTC) and Ultra-Reliable Low-Latency Communications
(URLLC). In this thesis, the coexistence of users with very low energy requirements, which are considered to belong to the mMTC category, and users with low latency requirements, which are
considered to belong to the URLLC category, is examined. Network slicing aims to meet their heterogeneous requirements, while the use of an edge server aims to achieve low communication delays. Meanwhile, both orthogonal (OMA) and non-orthogonal (NOMA) multiple access
schemes are examined. In particular, in the first case, a convex optimization problem is formulated, that aims to minimize the MEC offloading delay threshold for URLLC users, while
satisfying the heterogeneous requirements of both sets of devices. Then, the optimal pairing and resource allocation for the considered traffic categories are investigated, exploiting the NOMA
scheme. URLLC devices are orthogonally assigned to subchannels, but they share their resource blocks with mMTC devices and an optimization problem is formulated, that aims to minimize
the energy consumption for the data offloading of mMTC users while satisfying the diverse requirements of both sets of devices. The original problem is decomposed into two sub-problems
and these sub-problems are solved iteratively so as to obtain an optimal solution. Finally, the effectiveness of the proposed methods is verified through simulations.

In particular, for the first part, the convex optimization problem to be solved is:
<br />
![ScreenShot](https://github.com/VasilikiZarkadoula/Diploma_Thesis-Network_Slicing_in_MEC_systems/blob/main/P1.PNG)
<br />
, where u denotes the URRLC device, m the MMTC device and k both types of devices.
<br />
For the second part, the optimization problem to be solved is:
<br />
![ScreenShot](https://github.com/VasilikiZarkadoula/Diploma_Thesis-Network_Slicing_in_MEC_systems/blob/main/P2.PNG)
<br />
which is decomposed into the following sub-problems:
<br />
![ScreenShot](https://github.com/VasilikiZarkadoula/Diploma_Thesis-Network_Slicing_in_MEC_systems/blob/main/P3.PNG)
![ScreenShot](https://github.com/VasilikiZarkadoula/Diploma_Thesis-Network_Slicing_in_MEC_systems/blob/main/P4.PNG)
<br />

<br />
The simulation results are also provided in the 'Graphs' folders.
