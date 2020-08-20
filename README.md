# StringENT

StringENT is an extension of the randomness test suite [ENT](https://www.fourmilab.ch/random/), which assesses the randomness of PRNGs and TRNGs aimed for cryptography uses. 

Our extension provides p-values for tests already included in the ENT suite, thus helping perform hypothesis testing with these, while also adding other tests which provide further robustness, and detect further flaws of randomness. Moreover, the StringENT is nearly as fast and efficient as ENT, that is, many factors faster than the official implementations of NIST SP 800-22 and Dieharder.

## Usage
Use the extension as if using [ENT](https://www.fourmilab.ch/random/). It is the output that is changed. P-values which are significantly close to zero convey low randomness. Remember that many sequences of an RNG need to be tested in order to assess that RNG. 

## Tests included 
So far the suite has the following tests:
-Entropy (already present in ENT, not too useful for critical use).
-Chi square test of equidistribution (already in ENT).
-Arithmetic mean test (already in ENT, but p-values added).
-Monte Carlo pi estimation test (already in ENT, but p-values added).
-Serial Correlation test (already in ENT, but p-values added).
-Runs test (new, but might be removed in the future due to correlations).
-Local means test (new, generalization of NIST's block frequency test).

In StringENT, all tests bar the Entropy test provide a p-value with which to try and reject the null hypothesis at confidence level $\alpha$. The null hypothesis is that the sequence provides from a discrete uniform distribution $\mathcal{U}\{0,255\}$.


## Useful bibliography
-[NIST SP 800-22](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf)
-[The Art of Computer Programming, Volume 2](https://dl.acm.org/doi/book/10.5555/270146)

## Authors
This work was done by Bittor Alaña as part of his final year thesis in Mathematics at UCM. It has been supervised by Elena Almaraz and Javier García from UCM, and partly by Julio Hernández-Castro from University of Kent and Darren Hurley-Smith, from Royal Holloway University.