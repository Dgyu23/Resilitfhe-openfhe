# ResiliTFHE

This repository contains the research prototype for the paper:

> **ResiliTFHE: A Resilient Framework for Torus Fully Homomorphic Encryption Evaluation on CPUs**

ResiliTFHE is a research prototype implemented based on OpenFHE (https://github.com/openfheorg/openfhe-development) v1.4.2. It extends the OpenFHE codebase with fault-injection support and operator-level error-detection and protection mechanisms, together with the evaluation code used to study them.


## Git Branches

- `main` corresponds to the OpenFHE base code together with the shared ResiliTFHE documentation. It does not contain evaluation-specific code.
- `evaluation` contains the shared documentation and the evaluation implementations.

## Evaluation

All evaluations must be run from the `evaluation` branch. Run `git switch evaluation` before executing an experiment.

### PBS Operator Breakdown

Run the following commands to build and execute the PBS operator-breakdown evaluation:

```bash
cmake -S . -B build-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_EXAMPLES=ON \
  -DWITH_PBS_OPERATOR_BREAKDOWN=ON

cmake --build build-release --target evaluation-pbs-operator-breakdown -j

chmod +x scripts/evaluation_pbs_operator_breakdown.sh

scripts/evaluation_pbs_operator_breakdown.sh \
  build-release/bin/examples/binfhe/evaluation-pbs-operator-breakdown \
  results/pbs-operator-breakdown
```

This evaluation runs PBS across the predefined parameter configurations and records the operator-level timing breakdown. Inspect `results/pbs-operator-breakdown/summary.log` to verify that `NTT`, `INTT`, `VCMM`, `KeySwitch`, and `Decompose` collectively dominate the PBS execution time.

### PBS Operator Fault Injection

Run the following commands to build and execute the PBS operator fault-injection evaluation:

```bash
cmake -S . -B build-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_EXAMPLES=ON \
  -DWITH_PBS_SBF_FAULT_INJECTION=ON

cmake --build build-release --target evaluation-pbs-operator-fault-injection -j

chmod +x scripts/evaluation-pbs-operator-fault-injection.sh

scripts/evaluation-pbs-operator-fault-injection.sh \
  build-release/bin/examples/binfhe/evaluation-pbs-operator-fault-injection \
  results/pbs-operator-fault-injection
```

This evaluation injects transient single-bit-flip faults into the five PBS operators across the predefined parameter configurations and reports the decryption error rate (DER) in `results/pbs-operator-fault-injection/summary.log`. All five operators have non-negligible DER; `VCMM` and `NTT` are the most error-sensitive, with DER close to \(1 - 1/p\), where \(p\) denotes the plaintext modulus.
