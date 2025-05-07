# **Lecture: Introduction to Monte Carlo Methods**

## **1. Overview of Monte Carlo Methods**
Monte Carlo (MC) methods are computational algorithms that rely on repeated random sampling to obtain numerical results. They are particularly useful for:

- Solving problems with **high dimensionality**
- Evaluating **complex integrals**
- Simulating **stochastic processes**
- Performing **statistical physics** calculations

> *"Monte Carlo" refers to the famous casino, highlighting the role of randomness in these methods.*

---

## **2. Key Concepts**

### **A. Random Sampling**
- Generate random numbers from probability distributions
- Transform uniform random numbers to desired distributions using:
  - Inverse transform method
  - Rejection sampling
  - Markov Chain Monte Carlo (MCMC)

### **B. Law of Large Numbers**
- As sample size N → ∞, sample mean converges to expected value
- Basis for MC integration: ∫f(x)dx ≈ (1/N)Σf(xᵢ)

### **C. Central Limit Theorem**
- Provides error estimates for MC results
- Error decreases as 1/√N

---

## **3. Major Monte Carlo Techniques**

| Method | Description | Applications |
|--------|-------------|--------------|
| **Simple MC** | Direct random sampling | Integration, particle physics |
| **MCMC** | Construct Markov chain to sample from distribution | Bayesian statistics, statistical mechanics |
| **Metropolis-Hastings** | MCMC with acceptance criterion | Molecular simulations |
| **Gibbs Sampling** | Special case of MCMC for conditional distributions | Machine learning |
| **Importance Sampling** | Biased sampling for variance reduction | Rare event simulation |

---

## **4. Applications in Computational Science**

### **A. Statistical Physics**
- Ising model simulations
- Phase transition studies
- Polymer chain conformations

### **B. Financial Mathematics**
- Option pricing (Black-Scholes)
- Risk assessment
- Portfolio optimization

### **C. Computational Chemistry**
- Molecular dynamics (supplement to MD)
- Free energy calculations
- Protein folding

### **D. Machine Learning**
- Bayesian inference
- Training deep belief networks
- Reinforcement learning

---

## **5. Algorithm Example: Metropolis-Hastings**

```python
import numpy as np

def metropolis_hastings(target_dist, n_iters, initial_state):
    samples = [initial_state]
    current = initial_state
    
    for i in range(n_iters):
        proposal = current + np.random.normal()
        acceptance = min(1, target_dist(proposal)/target_dist(current))
        
        if np.random.rand() < acceptance:
            current = proposal
        samples.append(current)
    
    return samples
