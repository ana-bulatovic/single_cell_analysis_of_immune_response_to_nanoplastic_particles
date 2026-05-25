from shutil import which
import sys

if which("Rscript") is None:
    print("Rscript nije pronađen u PATH.")
    print("Instaliraj R (CRAN ili Posit) i proveri da li je Rscript dostupan u terminalu.")
    print("Primer: instaliraj R, zatvori terminal, pa pokreni `where Rscript`.")
    print("Ako se pojavi putanja, onda možeš koristiti `Rscript scripts/azimuth_annotation.R`.")
    sys.exit(1)

print("Rscript je pronađen. Možeš pokrenuti: Rscript scripts/azimuth_annotation.R")
