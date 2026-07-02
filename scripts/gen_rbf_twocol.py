"""
gen_rbf_twocol.py — Genera la versión doble columna (revtex4-2, estilo
APS/PRE) a partir del borrador principal, para mantener ambas versiones
sincronizadas.

Lee  docs/borrador_rebotes_gravedad.tex
Escribe docs/borrador_rbf.tex

Transformaciones:
  - preámbulo revtex4-2 con \\author/\\affiliation/\\email
  - abstract y \\keywords antes de \\maketitle (convención revtex)
  - figure → figure* y table → table* (todas las figuras/tablas son más
    anchas que una columna), [htbp] → [t]
  - sección Agradecimientos → entorno acknowledgments

Uso:  python3 scripts/gen_rbf_twocol.py && (cd docs && pdflatex borrador_rbf.tex)
"""

import re
from pathlib import Path

SRC = Path("docs/borrador_rebotes_gravedad.tex")
DST = Path("docs/borrador_rbf.tex")

s = SRC.read_text()

# --- abstract y palabras clave --------------------------------------------------
m = re.search(r"\\begin{abstract}(.*?)\\end{abstract}", s, re.S)
abstract_full = m.group(1)
abstract_body, _, tail = abstract_full.partition("\\medskip")
mkw = re.search(r"Palabras clave:}\s*(.*?)\s*$", tail.strip(), re.S)
keywords = re.sub(r"\s+", " ", mkw.group(1)).rstrip(".") if mkw else ""
# quitar las líneas de comentario decorativas
abstract_body = "\n".join(l for l in abstract_body.splitlines()
                          if not l.strip().startswith("% ==="))

# --- cuerpo: desde Introducción hasta la bibliografía ---------------------------
i0 = s.index(r"\section{Introducción}")
i1 = s.index(r"\end{thebibliography}") + len(r"\end{thebibliography}")
body = s[i0:i1]

# Agradecimientos → entorno revtex
body = re.sub(
    r"\\section\*{Agradecimientos}\s*\n% =+\n(.*?)(?=\n% =+\n\\bibliographystyle)",
    r"\\begin{acknowledgments}\n\1\\end{acknowledgments}\n",
    body, flags=re.S)
body = body.replace("\\bibliographystyle{unsrt}\n", "")

# figuras y tablas a ancho completo (dos columnas)
body = re.sub(r"\\begin{figure}\[\w+\]", r"\\begin{figure*}[tp]", body)
body = body.replace("\\end{figure}", "\\end{figure*}")
body = re.sub(r"\\begin{table}\[\w+\]", r"\\begin{table*}[tp]", body)
body = body.replace("\\end{table}", "\\end{table*}")

# figuras con rótulos en inglés (convención de la revista)
body = re.sub(r"(fig_[A-Za-z0-9_]+?)(_en)?\.(pdf|png)", r"\1_en.\3", body)

# --- documento revtex ------------------------------------------------------------
DOC = r"""% =============================================================================
% VERSIÓN RBF/BJP – Doble columna (revtex4-2, estilo APS/PRE)
% GENERADO por scripts/gen_rbf_twocol.py a partir de
% borrador_rebotes_gravedad.tex — NO editar a mano; editar el borrador
% principal y regenerar.
% Compilar: pdflatex borrador_rbf.tex  (dos pasadas)
% =============================================================================

\documentclass[aps,pre,twocolumn,reprint,superscriptaddress]{revtex4-2}

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[spanish,es-nodecimaldot]{babel}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{hyperref}

% --- Macros físicas ----------------------------------------------------------
\newcommand{\cor}{e}
\newcommand{\coR}[1]{e_{#1}}
\newcommand{\Tstar}{T^*}
\newcommand{\eps}{\varepsilon}
\newcommand{\epsps}{\varepsilon_{ps}}
\newcommand{\gstar}{g^*}
\newcommand{\lammps}{\textsc{Lammps}}
\newcommand{\TODO}[1]{\textbf{[TODO: #1]}}

% revtex ignora babel para este título de sección
\makeatletter
\AtBeginDocument{\@ifundefined{acknowledgmentsname}{}%
  {\renewcommand{\acknowledgmentsname}{Agradecimientos}}}
\makeatother

\begin{document}

\title{Partición de energía en el impacto de nanocúmulos Lennard-Jones
sobre superficies FCC: cuatro canales, escalado con el tamaño del
proyectil y adhesión}

% alternativa de correo: felipe.wasaff@uchile.cl
\author{Felipe Wasaff}
\email{felipe.wasaff.g@mail.pucv.cl}
\affiliation{Instituto de Matemáticas, Pontificia Universidad Católica de
Valparaíso, Valparaíso, Chile}

\author{Gonzalo Gutiérrez}
\affiliation{Departamento de Física, Facultad de Ciencias,
Universidad de Chile, Santiago, Chile}

\date{\today}

\begin{abstract}
__ABSTRACT__
\end{abstract}

\keywords{__KEYWORDS__}

\maketitle

__BODY__

\end{document}
"""

out = (DOC.replace("__ABSTRACT__", abstract_body.strip())
          .replace("__KEYWORDS__", keywords)
          .replace("__BODY__", body.strip()))
DST.write_text(out)
print(f"Generado: {DST} ({len(out.splitlines())} líneas)")
print(f"  keywords: {keywords}")
