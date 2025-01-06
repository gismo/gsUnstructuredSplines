# `gsUnstructuredSplines`: Unstructured spline constructions for G+Smo

|CMake flags|```-DGISMO_OPTIONAL="<other submodules>;gsUnstructuredSplines"```|
|--:|---|
|License|[MPL 2.0](https://www.mozilla.org/en-US/MPL/2.0/)|
|OS support|Linux, Windows, macOS|
|Build Status|[![ci](https://github.com/gismo/gsUnstructuredSplines/actions/workflows/ci.yml/badge.svg)](https://github.com/gismo/gsUnstructuredSplines/actions/workflows/ci.yml)|
|Repository|[gismo/gismo/gsUnstructuredSplines](https://github.com/gismo/gsUnstructuredSplines)|
|Dependencies|[gismo/gismo](https://github.com/gismo/gismo)|
|Developer|[Pascal Weinmueller](https://github.com/weinmueller),[Hugo Verhelst](https://github.com/hverhelst),[Andrea Farahat](https://github.com/AndreaFarahat)|
|Maintainers|[pascal.weinmueller@mtu.de](mailto:pascal.weinmueller@mtu.de),[h.m.verhelst@tudelft.nl](mailto:h.m.verhelst@tudelft.nl)|
|Last checked|22-08-2024|

## Installation
```
cd path/to/build/dir
cmake . -DGISMO_OPTIONAL="<other submodules>;gsUnstructuredSplines"
make
```

---

## Module overview

The `gsUnstructuredSplines` module provides ready-to-use unstructured spline constructions for smooth multi-patch modelling. The module provides the following unstructured spline constructions:
- **Approximate $C^1$** (`gsApproxC1Spline`)
  > Weinmüller, P. (2022). Weak and approximate C1 smoothness over multi-patch domains in isogeometric analysis, [***PhD Thesis***](https://epub.jku.at/obvulihs/content/titleinfo/7811106)
  >
  > Weinmüller, P., & Takacs, T. (2022). An approximate C1 multi-patch space for isogeometric analysis with a comparison to Nitsche’s method. [***Computer Methods in Applied Mechanics and Engineering***, 401, 115592.](https://doi.org/10.1016/j.cma.2022.115592)
  >
  > Weinmüller, P., & Takacs, T. (2021). Construction of approximate $C^1$ bases for isogeometric analysis on two-patch domains. [***Computer Methods in Applied Mechanics and Engineering***, 385, 114017.](https://doi.org/10.1016/j.cma.2021.114017)

- **Analysis-Suitable $G^1$** (`gsC1SurfSpline`)
  > Farahat, A. (2023). Isogeometric Analysis with $C^1$-smooth functions over multi-patch surfaces, [***PhD Thesis***](https://epub.jku.at/obvulihs/id/8255939)
  >
  > Farahat, A., Verhelst, H. M., Kiendl, J., & Kapl, M. (2023). Isogeometric analysis for multi-patch structured Kirchhoff–Love shells. [***Computer Methods in Applied Mechanics and Engineering***, 411, 116060.](https://doi.org/10.1016/j.cma.2023.116060)
  >
  > Farahat, A., Jüttler, B., Kapl, M., & Takacs, T. (2023). Isogeometric analysis with C1-smooth functions over multi-patch surfaces. [***Computer Methods in Applied Mechanics and Engineering***, 403, 115706.](https://doi.org/10.1016/j.cma.2022.115706)

- **Almost - $C^1$** (`gsAlmostC1`)
<<<<<<< HEAD
- **Degenerate patches (D-Patches)** (`gsDPatch`)
- **Multi-Patch B-Splines with Enhanced Smoothness** (`gsMPBESSpline`)
  > Buchegger, F., Jüttler, B., & Mantzaflaris, A. (2016). Adaptively refined multi-patch B-splines with enhanced smoothness. [***Applied Mathematics and Computation***, 272, 159-172.](https://doi.org/10.1016/j.amc.2015.06.055)

## Implementation aspects
The general implementation of unstructured spline constructions is provided by the `gsMappedSpline` and `gsMappedBasis` classes. These classes define a global basis construction through a linear combination of local basis functions. The linear combination is stored in the `gsWeightMapper`. In general, a mapped basis is configured as follows:

**TO DO**

=======
  > Takacs, T. & Toshniwal, D. (2023). Almost-$C^1$ splines: Biquadratic splines on unstructured quadrilateral meshes and their application to fourth order problems. [***Computer Methods in Applied Mechanics and Engineering***, 403, 115640.](https://doi.org/10.1016/j.cma.2022.115640)

- **Degenerate patches (D-Patches)** (`gsDPatch`)
  > Toshniwal, D., Speleers, H. & Hughes, T. J. (2017). Smooth cubic spline spaces on unstructured quadrilateral meshes with particular emphasis on extraordinary points: Geometric design and isogeometric analysis considerations. [***Computer Methods in Applied Mechanics and Engineering***, 327, 411-458.](https://doi.org/10.1016/j.cma.2017.06.008)

- **Multi-Patch B-Splines with Enhanced Smoothness** (`gsMPBESSpline`)
  > Buchegger, F., Jüttler, B., & Mantzaflaris, A. (2016). Adaptively refined multi-patch B-splines with enhanced smoothness. [***Applied Mathematics and Computation***, 272, 159-172.](https://doi.org/10.1016/j.amc.2015.06.055)

- **Compairison of the first four methods**
  > Verhelst, H. M. and Weinmüller, P. and Mantzaflaris, A. and Takacs, T. and Toshniwal, D. (2024). A comparison of smooth basis constructions for isogeometric analysis. [***Computer Methods in Applied Mechanics and Engineering***](https://doi.org/10.1016/j.cma.2023.116659)

## Implementation aspects
The general implementation of unstructured spline constructions is provided by the `gsMappedSpline` and `gsMappedBasis` classes. These classes define a global basis construction through a linear combination of local basis functions. The linear combination is stored in the `gsWeightMapper`. In general, a mapped basis is configured as follows:

>>>>>>> main
## Examples

<details>
<summary>Biharmonic equation</summary>

For more information, see the (Doxygen page)[url] corresponding to this file

</details>

<details>
<summary>Kirchhoff-Love shell model</summary>

For more information, see the (Doxygen page)[url] corresponding to this file

</details>

## Contributing to this module

## Publications based on this module

### Journal articles
<<<<<<< HEAD
1. Verhelst, H. M., Weinmüller, P., Mantzaflaris, A., Takacs, T., & Toshniwal, D. (2023). A comparison of smooth basis constructions for isogeometric analysis. ***arXiv preprint arXiv:2309.04405***.
=======
1. Verhelst, H. M., Weinmüller, P., Mantzaflaris, A., Takacs, T., & Toshniwal, D. (2023). A comparison of smooth basis constructions for isogeometric analysis. [***Computer Methods in Applied Mechanics and Engineering***, 419, 116659.](https://doi.org/10.1016/j.cma.2023.116659)
>>>>>>> main
1. Farahat, A., Verhelst, H. M., Kiendl, J., & Kapl, M. (2023). Isogeometric analysis for multi-patch structured Kirchhoff–Love shells. [***Computer Methods in Applied Mechanics and Engineering***, 411, 116060.](https://doi.org/10.1016/j.cma.2023.116060)
1. Farahat, A., Jüttler, B., Kapl, M., & Takacs, T. (2023). Isogeometric analysis with C1-smooth functions over multi-patch surfaces. [***Computer Methods in Applied Mechanics and Engineering***, 403, 115706.](https://doi.org/10.1016/j.cma.2022.115706)
1. Weinmüller, P., & Takacs, T. (2022). An approximate C1 multi-patch space for isogeometric analysis with a comparison to Nitsche’s method. [***Computer Methods in Applied Mechanics and Engineering***, 401, 115592.](https://doi.org/10.1016/j.cma.2022.115592)
1. Weinmüller, P., & Takacs, T. (2021). Construction of approximate $C^1$ bases for isogeometric analysis on two-patch domains. [***Computer Methods in Applied Mechanics and Engineering***, 385, 114017.](https://doi.org/10.1016/j.cma.2021.114017)
1. Buchegger, F., Jüttler, B., & Mantzaflaris, A. (2016). Adaptively refined multi-patch B-splines with enhanced smoothness. [***Applied Mathematics and Computation***, 272, 159-172.](https://doi.org/10.1016/j.amc.2015.06.055)

### PhD Theses
<<<<<<< HEAD
1. Verhelst, H.M. (2024). Isogeometric analysis of wrinkling, [***PhD Thesis***]()
=======
1. Verhelst, H.M. (2024). Isogeometric analysis of wrinkling, [***PhD Thesis***](https://doi.org/10.4233/uuid:0e4c3644-31a4-4157-983d-bd001d91b8ca)
>>>>>>> main
1. Farahat, A. (2023). Isogeometric Analysis with $C^1$-smooth functions over multi-patch surfaces, [***PhD Thesis***](https://epub.jku.at/obvulihs/id/8255939)
1. Weinmüller, P. (2022). Weak and approximate C1 smoothness over multi-patch domains in isogeometric analysis, [***PhD Thesis***](https://epub.jku.at/obvulihs/content/titleinfo/7811106)
---

# Changelog

***

<<<<<<< HEAD
#### Geometries:

![plot](./readme/dictionary_geometries.png)
=======
### Geometries
>>>>>>> main
