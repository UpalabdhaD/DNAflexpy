# Citation

If DNAflexpy is useful in your work, please cite it.

## From the command line

The installed package prints its own citation, with the version filled in:

```bash
DNAflexpy --citation
```

That is the authoritative version — it reads the same metadata this page does,
so the two cannot drift apart.

## BibTeX

```bibtex
@software{DNAflexpy,
  author    = {Upalabdha Dey and Rajesh Yella and Aditya Kumar},
  title     = {DNAflexpy: DNA flexibility profiling from sequence},
  year      = {2026},
  publisher = {GitHub},
  version   = {0.3.0},
  url       = {https://github.com/UpalabdhaD/DNAflexpy}
}
```

Replace `version` with the release you actually used:

```python
import DNAflexpy
print(DNAflexpy.__version__)
```

## Plain text

> Upalabdha Dey, Rajesh Yella, Aditya Kumar (2026). *DNAflexpy: DNA
> flexibility profiling from sequence.*
> <https://github.com/UpalabdhaD/DNAflexpy>

## Please also cite the parameter tables

DNAflexpy is a tool for applying published experimental parameter tables; it is
not the source of the numbers. The thirteen tables in
[`lookupNEW.yaml`](features.md) come from the primary literature — DNase I
cutting frequency, nucleosome positioning preference, wedge angle, propeller
twist, stacking free energy and so on — and **whichever tables your analysis
used should be cited alongside this software**.

!!! note "Those references are not yet recorded in the repository"
    The parameter tables ship without their source citations attached. Adding
    them to `docs/features.md` — one reference per table — would let users cite
    correctly without hunting through the literature themselves, and would make
    the package considerably more citable. This is the single most valuable
    documentation contribution outstanding.

## License

BSD 3-Clause. See [LICENSE](https://github.com/UpalabdhaD/DNAflexpy/blob/main/LICENSE).
