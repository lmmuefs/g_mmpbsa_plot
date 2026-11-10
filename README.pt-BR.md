[🇺🇸 English](README.md) · **🇧🇷 Português**

# mmpbsa_plot

Gráficos de MM/PBSA para quem usa o **g_mmpbsa** (versão para o campo de força GROMOS).
Desenvolvido no **LMM (Laboratório de Modelagem Molecular), UEFS**.

O script varre as subpastas, detecta os arquivos de cada sistema e abre um menu para
escolher o que plotar. No fim grava um relatório `.txt` com as escolhas e os resultados.

## Requisitos

Python 3.8+ com `numpy`, `pandas` e `matplotlib`.

```bash
pip install numpy pandas matplotlib
```

## Estrutura esperada

Uma subpasta por sistema (réplica, ligante, etc.), com os arquivos de saída do `g_mmpbsa`:

```
meu_projeto/
├── sistema_A/
│   ├── binding_energy.xvg            # energia vs. tempo (tempo + vdW, elet., polar, apolar, total)
│   ├── energy_summary.csv            # médias ± desvio dos componentes
│   ├── residues_energy_summary.csv   # decomposição por resíduo (trajetória inteira)
│   ├── contrib_MM.dat                # por frame e por resíduo: vdW + eletrostática (opcional)
│   ├── contrib_pol.dat               # por frame e por resíduo: solvatação polar   (opcional)
│   └── contrib_apol.dat              # por frame e por resíduo: solvatação apolar  (opcional)
├── sistema_B/
│   └── ...
└── mmpbsa_plot.py
```

Não precisa ter os três arquivos em todas as pastas: o menu só oferece o que existe.
Energias são lidas em kJ/mol e mostradas em kcal/mol; o tempo, em ns.

## Uso

```bash
python mmpbsa_plot.py                 # menu interativo, pasta atual
python mmpbsa_plot.py /caminho/dados  # outra pasta raiz
```

O primeiro menu pergunta o **idioma** (Português do Brasil ou English US). A escolha vale
para tudo: menus e mensagens, textos dos gráficos, cabeçalhos das tabelas e relatório
`.txt`. Em português, os números dos eixos usam vírgula decimal. O último idioma usado
fica salvo no `plot_style.ini`. Para pular a pergunta: `--lang pt` ou `--lang en`.
(O texto do `--help` é sempre em inglês.)

Depois o menu pergunta, nesta ordem:

1. **Sistemas** a incluir (`1,3`, `1-3` ou `all`).
2. **Legenda** de cada sistema. O padrão é o nome da pasta; Enter mantém.
3. **Gráficos**: energia no tempo (total, vdW, eletrostática, polar, apolar),
   médias ± desvio, decomposição por resíduo.
4. **Faixa de tempo** em ns (sem padrão: você digita início e fim). Vale para os gráficos
   de energia no tempo **e** para a decomposição por resíduo (veja Observações).
5. **Resíduos**: Top 5, Top 10 ou Todos. *Todos não é recomendado*: com centenas de
   resíduos o gráfico de barras fica ilegível, então só se gera o perfil por posição
   na sequência e o CSV completo.
6. **Ligante**: o script detecta o ligante (último resíduo que não é aminoácido
   padrão) e pergunta se mantém no ranking. Mantido, ele costuma dominar o topo.
7. **Cores**: hex separados por vírgula (`#0072B2,#D55E00`) ou Enter para a paleta padrão.
8. **Organização** das pastas de saída (ver abaixo).
9. **Formatos** (png, pdf, svg), exportar CSV dos dados, e janela interativa.

### Saída

Por padrão, tudo vai para uma pasta nova `MM_PBSA_plots_AAAAMMDD_HHMMSS/`, criada a cada
execução (então nada é sobrescrito), organizada por tipo de gráfico:

```
MM_PBSA_plots_AAAAMMDD_HHMMSS/
├── energy_vs_time/     # energia em função do tempo (um gráfico por componente)
├── summary/            # médias ± desvio dos componentes
├── per_residue/        # decomposição por resíduo (um por sistema)
├── tables/             # CSVs: rankings, tabelas "média ± desvio", dados dos gráficos
└── mmpbsa_report_AAAAMMDD_HHMMSS.txt
```

Use `-o nome_da_pasta` para escolher a pasta de saída (se ela já existir, é reaproveitada,
e arquivos com o mesmo nome são sobrescritos).

Em português as pastas ficam `energia_vs_tempo/`, `medias/`, `por_residuo/` e `tabelas/`.
Outras organizações disponíveis: **por formato** (`png/`, `pdf/`, `svg/`, útil na hora
de submeter) ou **tudo na mesma pasta**.

### O que é o gráfico de médias (`energy_components_mean`)

Barras com a média de cada termo da energia de ligação, por sistema, lidas do
`energy_summary.csv`: vdW, eletrostática, solvatação polar, solvatação apolar e total
(total = soma dos quatro). A média é sobre **toda a trajetória** e a barra de erro é o
**desvio padrão entre frames** (a flutuação da energia), não a incerteza da média.
Isso precisa constar na legenda da figura. O termo entrópico (−TΔS) não está incluído.

## Editando legendas, rótulos e fonte

Na primeira execução é criado o `plot_style.ini` na pasta raiz. Edite e rode de novo:

```ini
[names]
# aceita mathtext do matplotlib (comentários só em linha própria)
sistema_A = $\it{Mtb}$PS - KS1

[text_pt]
energy_ylabel = ENERGIA (kcal/mol)
residue_title = Decomposição de energia por resíduo
comp_vdw = vDW

[font]
family = Arial
size = 12
```

Os textos ficam em dois blocos, `[text_en]` e `[text_pt]`, um por idioma.
Um `plot_style.ini` de versões anteriores é atualizado sozinho na primeira execução
(as legendas que você definiu são mantidas).
A fonte padrão é **Arial**. Se ela não estiver instalada (comum em Linux), o script usa
o Liberation Sans, que tem as mesmas métricas, e registra isso no relatório.

Para ajustes finos fora do script:

- **SVG e PDF** mantêm o texto editável. Abra no Inkscape ou Illustrator
  (no SVG a fonte precisa estar instalada na máquina de quem abre).
- **Janela interativa** (opção do menu): o botão *Figure options* da barra do
  matplotlib altera título, rótulos, limites e cores antes de salvar.
- **CSV** dos dados: para refazer o gráfico em Origin, Excel ou Grace.

## Uso sem menu

Qualquer opção pula a pergunta correspondente. Com `-y` nada é perguntado: ele usa os
padrões (todos os sistemas, todos os formatos, ligante removido), então `--plots` é
obrigatório, e `--time` e `--top` também quando esses gráficos forem escolhidos.

```bash
python mmpbsa_plot.py . --systems 1,3 --plots total,vdw,residues \
    --time 20 100 --top 10 --ligand drop --formats png,svg --csv \
    --lang pt --layout type -y
```

`python mmpbsa_plot.py --help` lista todas as opções.

## Observações

- O gráfico de médias usa o `energy_summary.csv`, ou seja, a **trajetória inteira**.
  As médias do relatório para a faixa de tempo escolhida são calculadas a partir do `.xvg`
  e podem diferir.
- **Decomposição por resíduo e faixa de tempo.** O `residues_energy_summary.csv` é uma
  média única sobre os frames usados na execução do `g_mmpbsa`, então não dá para cortá-lo
  por tempo. Quando a faixa escolhida é menor que a trajetória, o script recalcula média e
  desvio (n−1) de cada resíduo a partir dos arquivos por frame `contrib_MM.dat`,
  `contrib_pol.dat` e `contrib_apol.dat` (gerados pelo `g_mmpbsa` com `-mmcon`, `-pcon`,
  `-apcon`), usando só os frames dentro da faixa. O `contrib_MM.dat` já soma vdW +
  eletrostática, então nesse caso o gráfico e as tabelas mostram **4 componentes** (MM,
  Polar, Apolar, Total) em vez de 5. Se a faixa cobrir a trajetória inteira, usa-se o CSV
  (5 componentes). Se os `contrib_*.dat` não existirem, usa-se o CSV e aparece um aviso. O
  relatório informa a faixa usada nos resíduos. Para ter vdW e eletrostática separadas em
  uma sub-faixa, rode o `g_mmpbsa` de novo com `-b`/`-e` dessa faixa.
- A tabela `*_top_N_table.csv` traz os valores no formato `média ± desvio`.

## Autoria

- Késsia Souza Santos: LMM (Laboratório de Modelagem Molecular), UEFS
