\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{mybeamer}[CU Boulder Beamer Presentation Template by @polmes]
\LoadClass[aspectratio = 169]{beamer}
\usetheme{metropolis}

\metroset{sectionpage = progressbar}
\metroset{subsectionpage = progressbar}
\metroset{progressbar = frametitle}
\metroset{numbering = fraction}
\metroset{block = fill}

% Progress bar thickness
\setlength{\metropolis@titleseparator@linewidth}{2pt}
\setlength{\metropolis@progressonsectionpage@linewidth}{2pt}
\setlength{\metropolis@progressinheadfoot@linewidth}{2pt}

% CU colors
\definecolor{gold}{RGB}{207, 184, 124}
\definecolor{black}{RGB}{0, 0, 0}
\definecolor{dark}{RGB}{86, 90, 92}
\definecolor{light}{RGB}{162, 164, 163}
\definecolor{white}{RGB}{255, 255, 255}
\definecolor{myblue}{RGB}{	52, 152, 219}
\setbeamercolor{normal text}{fg = black, bg = white}
\setbeamercolor{footer}{fg = white, bg = black}
\setbeamercolor{alerted text}{fg = gold, bg = light}
\setbeamercolor{example text}{fg = myblue, bg = white}

% Presentation mode
\RequirePackage{pgfpages}
\setbeamertemplate{note page}[plain]
% \setbeameroption{show notes on second screen = right}

% Replace enumerate from \note[item] with itemize
\def\beamer@setupnote{%
    \gdef\beamer@notesactions{%
        \beamer@outsideframenote{%
            \beamer@atbeginnote%
            \beamer@notes%
            \ifx\beamer@noteitems\@empty\else
                \begin{itemize}\itemsep=0pt\parskip=0pt%
                    \beamer@noteitems%
                \end{itemize}%
            \fi%
            \beamer@atendnote%
        }%
        \gdef\beamer@notesactions{}%
    }
}

% TiZ
\usetikzlibrary{tikzmark}
\usetikzlibrary{shapes, arrows, calc}
\tikzstyle{block} = [rectangle, draw, text width = 40mm, text centered, rounded corners, minimum height = 1.5\baselineskip]
\tikzstyle{line} = [draw, -latex']  

% Position blocks anywhere
\RequirePackage[absolute, overlay]{textpos}
\setlength{\TPHorizModule}{\paperwidth}
\setlength{\TPVertModule}{\paperheight}

% Custom blocks
\RequirePackage{tcolorbox}
\tcbuselibrary{skins}
\newtcolorbox{cublock}[1][Title]{
    enhanced,
    colback = white,
    colbacktitle = myblue,
    coltitle = white,
    colframe = black!70,
    boxrule = 1pt,
    titlerule = 0pt,
    arc = 5pt,
    title = {\strut#1},
    fonttitle = \bfseries,
    fuzzy shadow = {0pt}{-2pt}{0pt}{0.4pt}{black}
}

% Units typesetting
\RequirePackage[binary-units=true]{siunitx}
\sisetup{per-mode=symbol} % fractioned units
\sisetup{range-phrase=\,--\,} % range dash
\sisetup{range-units=single} % range units

% Math
\RequirePackage{mathtools}
% \RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{nicefrac}
\RequirePackage{bm}
% \RequirePackage{commath}
\RequirePackage{url}
% \RequirePackage{cancel}
\RequirePackage{physics}
% \RequirePackage{breqn}
% \RequirePackage{mleftright}
\RequirePackage{dsfont}
% Tables
\RequirePackage{array}
\RequirePackage{tabularx}
\RequirePackage{adjustbox}
\RequirePackage{longtable}
\RequirePackage{ltablex}
\RequirePackage{tikz-cd} %para diagrama conmutatitvo
\RequirePackage{multicol}
\RequirePackage{multirow}
\RequirePackage{booktabs}
\RequirePackage{tcolorbox}
\RequirePackage{threeparttable}
\RequirePackage{colortbl}
\RequirePackage[all,cmtip]{xy}
\RequirePackage{lipsum} %Just in case


% Cites in footones
\RequirePackage[style=authortitle, bibstyle=authortitle,backend=bibtex,maxnames=1]{biblatex}
\setbeamerfont{footnote}{size=\tiny}
\newcommand{\mycite}[1]{\footnote[frame,1]{\cite{#1}}}

\newcolumntype{L}[1]{>{\raggedright\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{C}[1]{>{\centering\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{R}[1]{>{\raggedleft\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}

% Normal quotes
\RequirePackage{csquotes}
\MakeOuterQuote{"}

% More symbols
\RequirePackage{pifont}

% Captions
\RequirePackage{caption}
\RequirePackage{subcaption}
\captionsetup{labelformat = empty, skip = 10pt, justification = centering}
\captionsetup[subfigure]{labelformat = empty, skip = 10pt, justification = centering}

% Numbered ToC
\setbeamertemplate{section in toc}[sections numbered]
% \setbeamertemplate{subsection in toc}[subsections numbered]

% Extra slides as "appendix"
\RequirePackage{appendixnumberbeamer}
\renewcommand{\appendixname}{\texorpdfstring{\translate{Appendix}}{Appendix}}

% The opposite of '\uncover'
\def\beamer@startmycovered{%
  \def\opaqueness<##1>##2{%
    \only<##1>{%
      \beamer@actions{%
        \expandafter\xdef\csname beamer@oldcolorhook%
        \the\beamer@coveringdepth\endcsname{\beamer@colorhook}%
        \expandafter\xdef\csname beamer@oldpgfextension%
        \the\beamer@coveringdepth\endcsname{\beamer@pgfextension}%
        {\globalcolorstrue\colorlet{beamer@freeze\the\beamer@coveringdepth}{bg}}%
        \xdef\beamer@colorhook{!##2!beamer@freeze%
          \the\beamer@coveringdepth\beamer@colorhook}%
        \gdef\beamer@pgfextension{!##2opaque}%
        \color{.}%
      }%
      {%
        \xdef\beamer@colorhook{\csname beamer@oldcolorhook%
          \the\beamer@coveringdepth\endcsname}%
        \xdef\beamer@pgfextension{\csname beamer@oldpgfextension%
          \the\beamer@coveringdepth\endcsname}%
        \color{.}%
      }}}%
  \ifnum\beamer@slideinframe<\beamer@minimum%ok, at beginning
  {%
    \beamer@saveanother%
    \advance\beamer@minimum by-\beamer@slideinframe%
    \beamer@slideinframe=\beamer@minimum%
    \beamer@uncoverbeforeactions%
    \beamer@restoreanother%
  }%
  \else%
  {%
    \beamer@saveanother%
    \advance\beamer@slideinframe by-\beamer@minimum%
    \beamer@uncoverafteractions%
    \beamer@restoreanother%
  }%
  \fi%
  \beamer@do%
%  }%
}
\long\def\beamer@makemycovered#1{\beamer@startmycovered#1\beamer@endcovered}
\def\cover{%
\alt{\beamer@makemycovered}{\beamer@fakeinvisible}}
\def\c@slideinframe{\beamer@slideinframe}

% Footer
\setbeamertemplate{footline}{%
    \begin{beamercolorbox}[wd = \textwidth, ht = 0.5cm, sep = 0.02cm]{footer}
        \hfill
        \raisebox{0.22cm}{
            \insertshortauthor
            \quad
            \textbullet
            \quad
            \inserttitle
            \quad
            |
            \quad
            \insertsection
            \quad
            \textbullet
            \quad
            \usebeamertemplate{frame numbering}
        }
        \hfill
    \end{beamercolorbox}%
}

% Keep footer in appendix
\AtBeginDocument{%
    \apptocmd{\appendix}{%
        \setbeamertemplate{footline}{%
        \begin{beamercolorbox}[wd = \textwidth, ht = 0.5cm, sep = 0.02cm]{footer}
            \hfill
            \raisebox{0.22cm}{
                \insertshortauthor
                \quad
                \textbullet
                \quad
                \inserttitle
                \quad
            }
            \hfill
            % \includegraphics[height = 0.6cm]{img/aems.pdf}
        \end{beamercolorbox}%
    }
    }{}{}
}


% Add footer to section pages
\AtBeginSection{
    {
        \setbeamertemplate{footline}{
           \begin{beamercolorbox}[wd = \textwidth, ht = 0.5cm, sep = 0.02cm]{footer}
                \hfill
                \raisebox{0.22cm}{
                    \insertshortauthor
                    \quad
                    \textbullet
                    \quad
                    \inserttitle
                }
                \hfill
                % \includegraphics[height = 0.6cm]{img/aems.pdf}
            \end{beamercolorbox}%
        }
        \frame{
            \sectionpage
        }
        \addtocounter{framenumber}{-1}
    }
}

% Add footer to subsection pages
\AtBeginSubsection{
    {
        \setbeamertemplate{footline}{
           \begin{beamercolorbox}[wd = \textwidth, ht = 1cm, sep = 0.2cm]{footer}
                \hfill
                \raisebox{0.22cm}{
                    \insertshortauthor
                    \quad
                    \textbullet
                    \quad
                    \inserttitle
                }
                \hfill
            \end{beamercolorbox}%
        }
        \frame{
            \sectionpage
        }
        \addtocounter{framenumber}{-1}
    }
}

% If one wants to modify the title page... here's the default code
\setbeamertemplate{title page}{%
    \begin{minipage}[b][\paperheight]{\textwidth}%
        \ifx\inserttitlegraphic\@empty\else\usebeamertemplate*{title graphic}\fi%
        \vfill%
        \ifx\inserttitle\@empty\else\usebeamertemplate*{title}\fi%
        \ifx\insertsubtitle\@empty\else\usebeamertemplate*{subtitle}\fi%
        \usebeamertemplate*{title separator}%
        \ifx\beamer@shortauthor\@empty\else\usebeamertemplate*{author}\fi%
        \ifx\insertdate\@empty\else\usebeamertemplate*{date}\fi%
        \ifx\insertinstitute\@empty\else\usebeamertemplate*{institute}\fi%
        \vfill%
        \vspace*{1mm}%
    \end{minipage}%
}

% Add footer to title page
\newcommand*{\copyrightnotice}[1]{\def\@copyrightnotice{#1}}
\renewcommand{\maketitle}[1][\item Welcome, everyone!]{
    \begingroup
        \setbeamercolor{normal text}{fg = white, bg = black}
        \setbeamertemplate{footline}{
            \begin{beamercolorbox}[wd = \textwidth, ht = 1.2cm, sep = 0.4cm]{footer}
                \hfill
                \ifdefined\@copyrightnotice
                    \raisebox{0.22cm}{
                        \begin{minipage}{0.65\paperwidth}
                            \begin{flushright}
                                Copyright~\copyright~by \insertshortauthor, \insertshortinstitute
                                \ifx\@copyrightnotice\@empty
                                \else
                                    \\
                                    \@copyrightnotice
                                \fi
                            \end{flushright}
                        \end{minipage}
                    }
                \fi
                % \includegraphics[height = 0.8cm]{img/aems.pdf}
            \end{beamercolorbox}
        }
        \begin{frame}
            \usebeamertemplate{title page}
            \vspace*{-50pt}
        \end{frame}
        \note[itemize]{#1}
        \addtocounter{framenumber}{-1}
    \endgroup
}
\RequirePackage{setspace}
% Table of Contents with "Outline" section name
\newcommand{\cutoc}{
    {
        \metroset{sectionpage = none}
        \section*{Estructura de la presentación}
    
        \begin{frame}{Estructura de la presentación}
            \renewcommand{\baselinestretch}{0.75}\normalsize
            \tableofcontents[hideallsubsections]
            \renewcommand{\baselinestretch}{1.0}\normalsize
        \end{frame}
    }
}

% Avoid overfull \vbox warnings
\vfuzz = 20pt


%This makes \pause work inside align, gather, and such
\let\save@measuring@true\measuring@true
\def\measuring@true{%
  \save@measuring@true
  \def\beamer@sortzero##1{\beamer@ifnextcharospec{\beamer@sortzeroread{##1}}{}}%
  \def\beamer@sortzeroread##1<##2>{}%
  \def\beamer@finalnospec{}%
}

%==============================
%========== Comandos ==========
%==============================
\newcommand{\mcU}{\mathcal{U}}
\newcommand{\mcV}{\mathcal{V}}
\newcommand{\mcO}{\mathcal{O}}
\newcommand{\mcI}{\mathcal{I}}
\newcommand{\mcL}{\mathcal{L}}
\newcommand{\mcS}{\mathcal{S}}
\newcommand{\hilbert}{{\sf H}}
\newcommand{\mcB}{\mathcal{B}}
\newcommand{\mcH}{\mathcal{H}}
\newcommand{\mcF}{\mathcal{F}}
\newcommand{\mcC}{\mathcal{C}}
\newcommand{\mcT}{\mathcal{T}}
\newcommand{\mcE}{\ensuremath{\mathcal{E}} }
\newcommand{\mcG}{\ensuremath{\mathcal{G}} }
\newcommand{\mcM}{\mathcal{M}}
\newcommand{\mcN}{\mathcal{N}}
\newcommand{\nnn}{\mathcal{N}}
\newcommand{\mmm}{\mathcal{M}}
\newcommand{\sss}{\mathcal{S}}
\newcommand{\mcD}{\mathcal{D}}
\newcommand{\mcA}{\mathcal{A}}
\newcommand{\mcP}{\mathcal{P}}
\newcommand{\rmi}{\text{i}}
\newcommand{\ie}{i.e.}
\newcommand{\avg}{\text{avg}}
\newcommand{\ef}{\text{ef}}
\newcommand{\Complex}{\mathbb{C}} %Para escribir al espacio de hilbert complejo
\newcommand{\Id}{\mathds{1}}% Para escribir el op. indentidad con notación chida
\newcommand{\CG}[1]{\mcC\left[#1\right]}
\newcommand{\Fuzzy}[1]{\mcF\left[#1\right]}
\newcommand{\pauli}[1]{\sigma_{#1}} %Para las matrices de pauli
\newcommand{\paulivec}[1]{\hat{#1}\cdot\vec{\sigma}} %Para vectores de Pauli unitarios
\newcommand{\cnot}{\text{C}_{\text{X}}} %Para el CNOT
\newcommand{\purity}[1]{\text{Pu}(#1)} %la pureza
\newcommand{\Abs}{\text{abs}} %abs
\newcommand{\rfroml}{f} %la función r(lambda)
\newcommand{\densityspace}[1]{\mcS(\hilbert_{#1})} %el espacio de operadores de densidad
\newcommand{\unitaryspace}[1]{\text{U}(#1)} %el espacio de operadores unitarios
\newcommand{\obspace}[1]{\mcL(\hilbert_{#1})} %el espacio de observables
\newcommand{\avgass}{\mcA_{\mcC}^{\avg}}
\newcommand{\maxass}{\mcA_{\mcC}^{\max}}
\DeclareMathOperator*{\Motimes}{\text{\raisebox{0.25ex}{\scalebox{0.8}{$\bigotimes$}}}} %para los productos tensoriales

\newcommand{\BlochSphere}{            \begin{tikzpicture}[line cap=round, line join=round, >=Triangle]
    \clip(-2.19,-2.49) rectangle (2.66,2.58);
    \draw [shift={(0,0)}, lightgray, fill, fill opacity=0.1] (0,0) -- (56.7:0.4) arc (56.7:90.:0.4) -- cycle;
    \draw [shift={(0,0)}, lightgray, fill, fill opacity=0.1] (0,0) -- (-135.7:0.4) arc (-135.7:-33.2:0.4) -- cycle;
    \draw(0,0) circle (2cm);
    \draw [rotate around={0.:(0.,0.)},dash pattern=on 3pt off 3pt] (0,0) ellipse (2cm and 0.9cm);
    \draw (0,0)-- (0.70,1.07);
    \draw [->] (0,0) -- (0,2);
    \draw [->] (0,0) -- (-0.81,-0.79);
    \draw [->] (0,0) -- (2,0);
    \draw [dotted] (0.7,1)-- (0.7,-0.46);
    \draw [dotted] (0,0)-- (0.7,-0.46);
    \draw (-0.08,-0.3) node[anchor=north west] {$\varphi$};
    \draw (0.01,0.9) node[anchor=north west] {$\theta$};
    \draw (-1.01,-0.72) node[anchor=north west] {$\hat{x}$};
    \draw (2.07,0.3) node[anchor=north west] {$\hat{y}$};
    \draw (-0.5,2.6) node[anchor=north west] {$\hat{z}=\ket{0}$};
    \draw (-0.4,-2) node[anchor=north west] {$-\hat{z}=\ket{1}$};
    \draw (0.4,1.65) node[anchor=north west] {$\ket{\psi}$};
    \scriptsize
    \draw [fill] (0,0) circle (1.5pt);
    \draw [fill] (0.7,1.1) circle (0.5pt);
  \end{tikzpicture}}

  \newcommand{\BlochSphereDensity}{            \begin{tikzpicture}[line cap=round, line join=round, >=Triangle]
    \clip(-2.19,-2.49) rectangle (2.66,2.58);
    \draw(0,0) circle (2cm);
    \draw [rotate around={0.:(0.,0.)},dash pattern=on 3pt off 3pt] (0,0) ellipse (2cm and 0.9cm);
    \draw [->] (0,0) -- (-0.81,-0.79); %Eje x
    \draw [->] (0,0) -- (2,0); %Eje y
    \draw [->] (0,0) -- (0,2); %Eje z
    \draw (0,0) -- (0.70,1.07); %Vector rho
    \draw [dotted] (0.0,1.46)-- (0.7,1.07); %Proyección en z
    \draw [dotted] (0,0)-- (0.7,-0.46); % Proyección en xy
    \draw [dotted] (0.7,1.07)-- (0.7,-0.46); %Proyección en xy
    \draw [dotted] (0.7,-0.46)-- (1.18,0.0); %Proyección en y
    \draw [dotted] (-0.47,-0.46)-- (0.7,-0.46); %Proyección en x
    \draw (-1.1,-0.15) node[anchor=north west] {$r_{x}$};
    \draw (1.,0.4) node[anchor=north west] {$r_{y}$};
    \draw (-0.5,1.75) node[anchor=north west] {$r_{z}$};
    \draw (-1.01,-0.72) node[anchor=north west] {$\hat{x}$};
    \draw (2.07,0.3) node[anchor=north west] {$\hat{y}$};
    \draw (-0.5,2.6) node[anchor=north west] {$\ket{0}$};
    \draw (-0.4,-2) node[anchor=north west] {$\ket{1}$};
    \draw (0.4,1.65) node[anchor=north west] {$\vec{r}_{\rho}$};
    \scriptsize
    \draw [fill] (0,0) circle (1.5pt);
    \draw [fill] (0.7,1.07) circle (0.5pt);
  \end{tikzpicture}}