# An Adaptive (CFAR) Detection Algorithm 

This repository focuses on the implementation and analysis of an **adaptive detection algorithm** as presented in the seminal paper [An Adaptive Detection Algorithm](https://ieeexplore.ieee.org/document/4104190). The repository is a work in progress, and we are continuously adding more simulations and analyses.

## Key Topics

### 1. CFAR Detectors
- **Problem Addressed**:
  - Detection of signals in the presence of unknown Gaussian noise using statistical hypothesis testing.
  - Signal detection is performed in one data vector, with another independent set of signal-free data vectors available. These vectors share the unknown covariance matrix of the noise in the former vector.
- **Solution**:
  - A likelihood ratio decision rule (LRT) is derived to address this problem.
  - The performance of the LRT is evaluated under:
    - Noise-only scenarios
    - Signal-plus-noise scenarios

### 2. Comparison of Detection Techniques
- **Simulations**:
  - A comprehensive comparison is made between the **Likelihood Ratio Test (LRT)** and the **Matched Filter**.
- **Key Insights**:
  - Performance metrics and trade-offs between the two approaches are analyzed to highlight their strengths and weaknesses.

---

We are actively expanding this repository with additional simulations, analyses, and implementations. Contributions and feedback are welcome!
