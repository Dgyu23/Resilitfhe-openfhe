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

### PBS Operator Protection

Run the following commands to build and execute the PBS operator-protection evaluation:

```bash
cmake -S . -B build-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_EXAMPLES=ON

cmake --build build-release --target evaluation-operator-protection -j

chmod +x scripts/evaluation-operator-protection.sh

scripts/evaluation-operator-protection.sh build-release/bin/examples/binfhe/evaluation-operator-protection results/operator-protection ntt
scripts/evaluation-operator-protection.sh build-release/bin/examples/binfhe/evaluation-operator-protection results/operator-protection intt
scripts/evaluation-operator-protection.sh build-release/bin/examples/binfhe/evaluation-operator-protection results/operator-protection vcmm
scripts/evaluation-operator-protection.sh build-release/bin/examples/binfhe/evaluation-operator-protection results/operator-protection decompose
scripts/evaluation-operator-protection.sh build-release/bin/examples/binfhe/evaluation-operator-protection results/operator-protection keyswitch
```

This evaluation details the execution latency and protection overhead of the five performance-critical PBS operators—NTT, INTT, VCMM, Decompose, and KeySwitch—which collectively account for most PBS execution time, across four schemes: the unprotected baseline, DMR, Naive ABFT, and ResiliTFHE. It validates our scheme's protection-overhead optimization relative to DMR and Naive ABFT.

### Protected Operator Fault-Detection Reliability

Run the following commands to build and execute the operator-level fault-detection evaluation:

```bash
cmake -S . -B build-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_EXAMPLES=ON

cmake --build build-release --target evaluation-operator-reliability -j

chmod +x scripts/evaluation-operator-reliability.sh

scripts/evaluation-operator-reliability.sh build-release/bin/examples/binfhe/evaluation-operator-reliability results/operator-reliability ntt
scripts/evaluation-operator-reliability.sh build-release/bin/examples/binfhe/evaluation-operator-reliability results/operator-reliability intt
scripts/evaluation-operator-reliability.sh build-release/bin/examples/binfhe/evaluation-operator-reliability results/operator-reliability vcmm
scripts/evaluation-operator-reliability.sh build-release/bin/examples/binfhe/evaluation-operator-reliability results/operator-reliability decompose
scripts/evaluation-operator-reliability.sh build-release/bin/examples/binfhe/evaluation-operator-reliability results/operator-reliability keyswitch
```

This evaluation injects transient single-bit-flip faults into the protected NTT, INTT, VCMM, Decompose, and KeySwitch implementations. For each selected operator, the script automatically evaluates every current fault site, STD128/STD192/STD256, and both single-fault (1F) and two-fault (2F) injections; in 2F experiments, the two faults target distinct dynamic occurrences of the same site during one invocation. The resulting `results/operator-reliability/<OPERATOR>/summary.csv` reports detected, undetected, and masked samples, together with detection and masking rates. A masked sample has no alarm and exactly matches the fault-free output, so it records a fault that did not affect the operator result rather than an undetected erroneous output.
