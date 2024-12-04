# NAQRA

_**N**eon **A**ccelerated **QR** **A**lgorithm_

## Introduction

An implementation of the [_QR Algorithm_](https://en.wikipedia.org/wiki/QR_algorithm) with _Double Wilkinson's Shift_ utilizing [_Neon_](https://developer.arm.com/Architectures/Neon) intrinsics, written in `C23`.

## Table of Contents

- [Introduction](#introduction)
- [Table of Contents](#table-of-contents)
- [Setup](#setup)
    - [Cloning the Repository](#cloning-the-repository)
    - [Compiling the Tests](#compiling-the-tests)
- [Usage](#usage)
    - [Interface](#interface)

## Setup

### Cloning the Repository

Clone the repository from [here](https://github.com/diantonioandrea/NAQRA):

```bash
git clone git@github.com:diantonioandrea/NAQRA.git
```

### Compiling the Tests

Compile the tests located under `src/`:

```bash
make
```

The executables are then located under `executables/`.

Tests require the size of the test matrix, along with the range of random values used to initialize it, specified as two real numbers. An optional argument can be provided to enable the output.

## Usage

Every method developed in **NAQRA** follows a structured naming convention with three parts, separated by underscores:

1. **Function**: Describes the primary operation or behavior of the method.
2. **Input(s)**: Specifies the type or nature of the inputs.
3. **Output**: Indicates the type or nature of the resulting output.

All methods are documented directly within the code for clarity and ease of use.

### Interface[^Interface]

Moreover, the repository provides an interface that includes the `Vector` and `Matrix` structures, some output methods, and the `Eigenvalues` function for higher-level usage. Vectors can be created, accessed, edited and deleted using the `NewVector`, `GetVectorAt`, `SetVectorAt`, and `FreeVector` methods. Similarly, matrices can be manipulated with methods whose name follows the same conventions.

[^Interface]: All tests are implemented with low level functions.