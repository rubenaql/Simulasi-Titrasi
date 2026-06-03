import os
import runpy

# Entry point deploy
# Pastikan path Simulasi.py di-load relatif terhadap folder ini.
_here = os.path.dirname(os.path.abspath(__file__))
runpy.run_path(os.path.join(_here, "Simulasi_v2.py"), run_name="__main__")

