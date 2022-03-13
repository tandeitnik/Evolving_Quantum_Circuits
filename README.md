#   Evolving Quantum Circuits

 <h4 align="center"> 
	🚧  READ ME UNDER CONSTRUCTION  🚧
</h4>

<!-- PROJECT LOGO -->
<br />
<div align="center">
  <img src="https://i.imgur.com/uCyUURi.png" alt="Logo">
  

  <p align="center">
    🧬 <b>Application of genetic algorithm to quantum computing</b> 🧬
    <br />
  </p>
</div>

Summary
=================
<!--ts-->
   * [About](#about)
   * [Publications](#publications)
   * [Requirements](#requirements)
   * [How-to use it](#how-to-use-it)
   * [List of functions](#list-of-functions)
      * [Stabilizer functions](#stabilizer-functions)
      * [Evolution functions](#evolution-functions)
      * [QECC functions](#qecc-functions)
      * [Utility functions](#utility-functions)
   * [Author](#author)

<!--te-->

## About

This repository stores the scripts made to explore genetic algorithms (GA) as a tool to evolve quantum circuits. GA mimic the evolutionary process of natural selection displayed by nature to solve optimization problems.

For a project for my master's degree in Physics, I raised the question of whether genetic algorithms have the potential to solve problems in the context of quantum computing, where human intuition decreases as physical systems grow. Specifically, we focused on the evolution of quantum error-correcting codes (QECCs) within the stabilizer code formalism. By specifying an appropriate fitness function, we showed that we could evolve celebrated codes, such as the Perfect, Color, and Shor's code with respectively 5, 7, and 9 qubits, in addition to new unanticipated examples. Additionally, we compared it with a brute force random search and verified the supremacy of GA over it. Given the results, we foresee that genetic algorithms can become valuable tools to perform complex applications in quantum systems and produce tailored circuits that satisfy restrictions imposed by hardware.

This readme brings a brief starter guide of utilizing the scripts contained in the repository. For a more in-depth explanation, please refer to the dissertation. Feel free to contact me via e-mail; I will gladly answer your questions. 

## Publications


🚧  Soon links to the dissertation and the article. 🚧 

## Requirements

 - [x] Python
 - [x] numpy
 - [x] qiskit
 - [x] itertools
 - [x] networkx
 - [x] basic knowledge of quantum computation, stabilizer formalism and quantum error correction codes

*This code is Python 2 compatible.*

## How-to use it



## List of functions

According to my classification, all functions are sorted into four different scripts for better organization. The following sub-sections are a comprehensive list of the functions with a brief description of their functionalities.


## Stabilizer functions

<!--ts-->
   * [randcirc](#randcirc)
   * [gene_express](#gene_express)  
   * [final_state](#final_state)  
<!--te-->

<hr>

### randcirc

```python
randcirc(N,T,adj_mat)
```
**Input:**  a) `N` , number of qubits; b) `T` , number of gates per qubit; c) `adj_mat` , adjacency matrix of the lattice.

**Output:** a) `circuit` , a random Clifford circuit array.

**Description:** This function generates a random Clifford circuit over `N` qubits. CNOT gates satisfies the neighbour structure imposed by the adjacency matrix `adj_mat`. In `T` loops, a random Clifford gate (I, H, P, CNOT) is selected for each qubit in the crescent order of indexes, thus the final circuit has NT gates.

<hr>

### gene_express

```python
gene_express(N,circ)
```
**Input:**  a) `N` , number of qubits; b) `circ` , Clifford circuit array.

**Output:** a) `circuit` , qiskit `QuantumCircuit` object.

**Description:** This function translates a Clifford circuit array into its qiskit equivalent `QuantumCircuit` object.

<hr>

### final_state

```python
final_state(N,circ)
```
**Input:**  a) `N` , number of qubits; b) `circ` , Clifford circuit array.

**Output:** a) `out_state` , qiskit result state vector.

**Description:** This function returns the output state of the computation of the Clifford circuit `circ` to the initial state |0>^N.

<hr>

## Evolution functions

## QECC functions

## Utility functions

## Author

 <img style="border-radius: 50%;" src="https://i.imgur.com/KuuiVV2.png" width="150px;" alt=""/>

**Daniel Ribas Tandeitnik**<br />
*Physicist and production engineer.* 

Feel free to get in touch with me.


[![Linkedin Badge](https://img.shields.io/badge/-Daniel_Tandeitnik-blue?style=flat-square&logo=Linkedin&logoColor=white&link=https://www.linkedin.com/in/tandeitnik/)](https://www.linkedin.com/in/tandeitnik/) [![Gmail Badge](https://img.shields.io/badge/-tandeitnik@gmail.com-c14438?style=flat-square&logo=Gmail&logoColor=white&link=mailto:tandeitnik@gmail.com)](mailto:tandeitnik@gmail.com)
