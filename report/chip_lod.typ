// chip_lod report. Every number, table and figure is read from results/ — recompile after a rerun.
// Compile from the repo root: typst compile --root . report/chip_lod.typ

#let res = "../results/python/"
#let s = json(res + "summary.json")
#let hc = json("../results/hash_check.json")
#let F = s.f
#let B = s.bayes_example

// ---------- helpers ----------

#let load(name) = {
  let rows = csv(res + name)
  let h = rows.first()
  rows.slice(1).map(r => h.zip(r).to-dict())
}
#let near(a, b) = calc.abs(float(a) - float(b)) <= 1e-9 * calc.max(1.0, calc.abs(float(b)))
#let isna(x) = x == none or x == "NA"
#let fmt(x, d: 2) = if isna(x) { [--] } else { str(calc.round(float(x), digits: d)) }
#let pct(x, d: 2) = if isna(x) { [--] } else { str(calc.round(float(x) * 100, digits: d)) + "%" }
#let istr(x) = if isna(x) { [--] } else { str(calc.round(float(x))) }
#let sci(x, d: 1) = {
  let v = float(x)
  if v == 0 { return [0] }
  let e = calc.floor(calc.log(v, base: 10))
  let m = v / calc.pow(10.0, e)
  [#str(calc.round(m, digits: d)) × 10#super[#str(e)]]
}
#let prob(x) = {
  if isna(x) { return [--] }
  let v = float(x)
  if v == 0 { [0] } else if v == 1 { [1] } else if v < 0.001 { sci(v) } else if v > 0.9999 { [> 0.9999] } else { str(calc.round(v, digits: 4)) }
}
#let ppl(x) = { // a "number of people" in the natural-frequency examples
  let v = float(x)
  if v >= 10 { str(calc.round(v, digits: 0)) } else if v >= 0.01 { str(calc.round(v, digits: 1)) } else { sci(v) }
}
#let dfmt(x) = if x == none { [never (within the grid)] } else { [#x reads] }
#let keybox(body) = block(fill: luma(240), inset: 10pt, radius: 4pt, width: 100%, body)
#let note(body) = block(inset: (left: 8pt), stroke: (left: 2pt + luma(180)), body)

#set document(title: "chip_lod: how many reads make a CHIP call?", author: "chip_lod")
#set page(paper: "a4", margin: 2.2cm, numbering: "1")
#set text(size: 10.5pt)
#set par(justify: true)
#set heading(numbering: "1.")
#show table: set text(size: 9pt)
#show figure.caption: set text(size: 9pt)

// part titles (unnumbered)
#let part(title) = {
  pagebreak(weak: true)
  align(center, block(above: 6pt, below: 14pt, text(size: 16pt, weight: "bold", title)))
}

#align(center)[
  #text(size: 18pt, weight: "bold")[chip_lod: how many alternate reads make a CHIP call?] \
  #v(2pt)
  #text(size: 11pt)[A plain-language guide to the model, the analyses and how to read them]
]

#v(6pt)
#keybox[
  *The question.* At one known position (for example a _DNMT3A_ or _JAK2_ hotspot), a sequencer reports $n$ reads,
  and $k$ of them show the variant base. Some of those reads can be sequencing errors. How many variant reads do we
  need before we can trust the call, and how small a variant can we detect at all?

  *Part I (frequentist).* With #F.reference_depth reads and an error rate of #F.reference_error_rate, a call needs at least
  *#F.k_star variant reads*. With that rule a variant at *#pct(F.lod_vaf) VAF* or more is found in #pct(F.power, d: 0) of
  samples: that is the *limit of detection* (LoD). The LoD rises with the error rate (from #pct(F.by_error_rate.first().lod_vaf)
  at $e$ = #F.by_error_rate.first().error_rate to #pct(F.by_error_rate.last().lod_vaf) at $e$ = #F.by_error_rate.last().error_rate)
  and falls with depth. Detecting a #pct(F.required_depth.last().target_vaf, d: 0) variant needs
  #F.required_depth.last().min_depth reads, and a #pct(F.required_depth.first().target_vaf, d: 1) variant
  #F.required_depth.first().min_depth reads.

  *Part II (Bayesian).* The frequentist test says whether reads are surprising _if there were no variant_. The Bayesian
  analysis answers the question a clinician usually asks: _how likely is it that the variant is really there?_ It
  also accounts for how rare the variant is to begin with. At depth #s.a3.depth, believing in the variant with 95%
  certainty takes #s.a3.by_prior.at(0).k_h1 reads if variants are rare (1 in #calc.round(1 / s.a3.by_prior.at(0).prior)),
  and #s.a3.by_prior.at(1).k_h1 if they are a coin flip.
]

#outline(depth: 1, indent: 1em)

= The question

Clonal haematopoiesis of indeterminate potential (CHIP) is defined by a somatic mutation in blood at a variant
allele fraction (VAF) of at least #pct(s.model.p_thr, d: 0). Blood cells are diploid and the variant sits on one
of the two copies, so a VAF of 2% means about 4% of cells carry it.

Sequencing is not perfect. Even at a position where nobody carries the variant, a small fraction of reads show
the variant base because of sequencing or PCR errors. With an error rate $e$ = #F.reference_error_rate per base,
and a third of errors landing on any particular wrong base, a site with #F.reference_depth reads is expected to show
*#fmt(F.lambda_bg, d: 3) error reads* on average. So when we see 1, 2 or 5 variant reads, is that a variant or noise?

This report looks at the question in two ways. *Part I* uses only classical (frequentist) statistics: a
significance test against the error background, and the limit of detection that follows from it. *Part II* adds
the Bayesian view, explained from scratch.

All results come from the `chip_lod` package (Python, with an R port that produces byte-identical files). All
settings are in `config/params.json`.

// =====================================================================================================
#part[Part I — The frequentist analysis]

#let f1 = load("f1_lod_grid.csv")
#let f2 = load("f2_required_depth.csv")
#let f3 = load("f3_power_curve.csv")
#let f4 = load("f4_sim_power.csv")
#let f5 = load("f5_lod_function_example.csv")
#let f-errs = f1.map(r => r.error_rate).dedup()
#let f-depths = f1.map(r => r.depth).dedup()
#let ref-ex = f5.filter(r => int(r.n) == F.reference_depth)

= How the frequentist test works

The frequentist approach asks one question: *if there were no variant, how surprising would these reads be?*
It needs only a model of the errors. Here is the procedure at depth #F.reference_depth with $e$ = #F.reference_error_rate.

+ *Background.* Errors produce on average $lambda_0 = n dot e slash 3$ = #fmt(F.lambda_bg, d: 3) variant reads. The
  actual number varies from sample to sample and follows a Poisson distribution, the standard model for counts
  of rare, independent events.
+ *p-value.* The p-value of $k$ reads is the probability that errors alone produce $k$ _or more_ variant reads.
  #for r in ref-ex [With $k$ = #r.k the p-value is #prob(r.pvalue). ]
+ *Critical value $k^*$.* We call the variant when the p-value is at most $alpha$ = #F.alpha, i.e. when errors would
  produce that many reads in fewer than #calc.round(F.alpha * 100)% of variant-free samples. The smallest such
  count is $k^*$ = *#F.k_star*. Because read counts are whole numbers, the actual false-positive rate is
  #prob(F.actual_alpha), a little below $alpha$.
+ *Power.* A real variant at VAF $p$ produces variant reads at the rate $r(p) = p(1-e) + (1-p) e slash 3$ per read
  (correct reads from the mutated copy plus misread normal reads). The *power* is the chance that such a variant
  reaches $k^*$ reads and is called.
+ *Limit of detection.* The LoD is the smallest VAF that is called in at least #pct(F.power, d: 0) of samples.
  At depth #F.reference_depth this is *#pct(F.lod_vaf)*.

#note[
  *How to read this.* $k^*$ is a property of the depth and the error rate, and is fixed before looking at the data.
  The LoD is a property of the assay: it says which variants the assay reliably detects, not whether this
  particular sample carries one. A variant below the LoD can still be called, just not reliably; a sample
  without a call can still carry a variant below the LoD.
]

The rest of Part I varies the three ingredients one at a time (error rate, depth and VAF), focusing on VAFs up to
#pct(F.vaf_max, d: 0), the range where calling is difficult. Throughout, $alpha$ = #F.alpha and the required power is
#pct(F.power, d: 0).

= Influence of the error rate

#figure(
  table(
    columns: 5,
    align: right,
    table.header([error rate $e$], [background $lambda_0$], [$k^*$], [actual false-positive rate], [LoD VAF]),
    ..F.by_error_rate.map(r => ([#r.error_rate], fmt(r.lambda_bg, d: 3), [#r.k_star], prob(r.actual_alpha),
      pct(r.lod_vaf))).flatten(),
  ),
  caption: [The effect of the error rate at depth #F.reference_depth (`f1_lod_grid.csv`).],
)

A higher error rate means more background reads, so a call needs more variant reads ($k^*$ rises from
#F.by_error_rate.first().k_star to #F.by_error_rate.last().k_star), and the LoD rises with it. Going from the
cleanest assay ($e$ = #F.by_error_rate.first().error_rate, typical of UMI/duplex consensus reads) to the noisiest
($e$ = #F.by_error_rate.last().error_rate) raises the LoD about
#calc.round(float(F.by_error_rate.last().lod_vaf) / float(F.by_error_rate.first().lod_vaf), digits: 1)-fold.

The "actual false-positive rate" column shows the effect of whole-number thresholds. The test promises at most
#pct(F.alpha, d: 0) false positives, but because $k^*$ jumps in whole steps, the real rate is anywhere between almost
0 and #pct(F.alpha, d: 0).

#figure(
  table(
    columns: 1 + f-errs.len(),
    align: right,
    table.header([depth], ..f-errs.map(e => [$e$ = #fmt(e, d: 4)])),
    ..f-depths.map(d => ([#d], ..f-errs.map(e => {
      let r = f1.find(r => r.depth == d and r.error_rate == e)
      [#pct(r.lod_vaf) (#r.k_star)]
    }))).flatten(),
  ),
  caption: [LoD VAF for every combination of depth and error rate, with $k^*$ in brackets (`f1_lod_grid.csv`).],
)

= Influence of the read depth

#figure(
  image("../results/figures/fig_f1_lod_depth.svg", width: 72%),
  caption: [Limit of detection against depth for six error rates (log–log). The dotted line is the 2% CHIP
    threshold.],
)

More reads make the variant reads stand out from the background, so the LoD falls with depth. While errors are
negligible ($k^*$ stays at 1 or 2), it falls almost in proportion to $1 slash n$: ten times the depth gives a ten
times smaller LoD. Once the background amounts to several reads, the variant has to rise above the random scatter
of those error reads, and the LoD falls more slowly, closer to $1 slash sqrt(n)$. At $e$ = #F.reference_error_rate the LoD goes from
#pct(f1.find(r => near(r.error_rate, F.reference_error_rate)).lod_vaf) at #f-depths.first() reads to
#pct(f1.filter(r => near(r.error_rate, F.reference_error_rate)).last().lod_vaf) at #f-depths.last() reads.
The lines are not perfectly smooth: each time the depth grows enough for $k^*$ to step up by one read, the LoD
briefly gets _worse_ before improving again. This is the "saw-tooth" effect of whole-number thresholds.

*How deep do I need to sequence?* The table below inverts the question. It gives the smallest depth at which a
variant of a given VAF is detected with #pct(F.power, d: 0) power.

#let f2-targets = f2.map(r => r.target_vaf).dedup()
#figure(
  table(
    columns: 1 + f2-targets.len(),
    align: right,
    table.header([error rate $e$], ..f2-targets.map(v => [VAF #pct(v, d: 1)])),
    ..f-errs.map(e => ([#fmt(e, d: 4)], ..f2-targets.map(v => {
      let r = f2.find(r => r.error_rate == e and r.target_vaf == v)
      [#istr(r.min_depth) ($k^*$ = #istr(r.k_star))]
    }))).flatten(),
  ),
  caption: [Smallest depth with #pct(F.power, d: 0) power to detect each VAF (`f2_required_depth.csv`).],
)

#figure(
  image("../results/figures/fig_f3_required_depth.svg", width: 66%),
  caption: [Depth needed for #pct(F.power, d: 0) power against error rate, for three target VAFs.],
)

#let nonmono = f2-targets.filter(v => {
  let ds = f-errs.map(e => f2.find(r => r.error_rate == e and r.target_vaf == v).min_depth).map(float)
  range(ds.len() - 1).any(i => ds.at(i + 1) < ds.at(i))
})
Small VAFs are expensive: halving the VAF roughly doubles (or more) the depth needed, and a noisy assay can
need several times more reads than a clean one.
#if nonmono.len() > 0 [
  One result looks paradoxical: for #nonmono.map(v => pct(v, d: 1)).join(", ", last: " and ") VAF, a
  slightly _higher_ error rate sometimes needs _fewer_ reads. This happens when both error rates use the same
  $k^*$. The extra error reads then also count towards the variant's total and help it cross the threshold.
  The price is a higher actual false-positive rate (previous section), so the noisier assay is not really better.
]

= Influence of the VAF

#figure(
  image("../results/figures/fig_f2_power.svg", width: 90%),
  caption: [Power (probability of a call) against true VAF up to #pct(F.vaf_max, d: 0). Left: depths at $e$ =
    #F.reference_error_rate. Right: error rates at depth #F.reference_depth. The dotted line is #pct(F.power, d: 0)
    power (`f3_power_curve.csv`).],
)

#let f3-ref = f3.filter(r => near(r.error_rate, F.reference_error_rate))
#let f3-depths = f3-ref.map(r => r.depth).dedup()
#let show-vafs = (0.001, 0.0025, 0.005, 0.01, 0.02).filter(v => v <= float(F.vaf_max) + 1e-12)
#figure(
  table(
    columns: 1 + f3-depths.len(),
    align: right,
    table.header([true VAF], ..f3-depths.map(d => [n = #d ($k^*$ = #f3-ref.find(r => r.depth == d).k_star)])),
    ..show-vafs.map(v => ([#pct(v, d: 2)], ..f3-depths.map(d => {
      let r = f3-ref.find(r => r.depth == d and near(r.true_vaf, v))
      if r == none [--] else [#pct(r.power, d: 1)]
    }))).flatten(),
  ),
  caption: [Probability of calling a variant of a given VAF at $e$ = #F.reference_error_rate (selected rows of
    `f3_power_curve.csv`).],
)

#let f1-ref = f1.filter(r => near(r.error_rate, F.reference_error_rate))
#let d01 = f1-ref.find(r => float(r.lod_vaf) <= 0.001)
For a fixed assay, power rises steeply with VAF: below the LoD a variant is called only some of the time, and above it
almost always. Between about 0.1% and 1% VAF, the region of interest for early clonal haematopoiesis, the answer
depends strongly on depth. At $e$ = #F.reference_error_rate, a variant at the
#pct(F.required_depth.last().target_vaf, d: 0) CHIP threshold is reliably detected from
#F.required_depth.last().min_depth reads onwards, #pct(F.required_depth.first().target_vaf, d: 1) needs
#F.required_depth.first().min_depth reads, and 0.1% #if d01 == none [is not reached at any depth in our grid
(up to #f-depths.last() reads)] else [is first reached at #d01.depth reads in our depth grid].

= A function for the frequentist LoD

`lod_frequentist(n, k, e, alpha, power)` applies the whole of Part I to a list of sites. Given vectors of depths
$n$ and variant-read counts $k$, it returns for every site: the background $lambda_0$, the p-value, $k^*$, whether the
site is called, and the LoD at that depth.

#grid(columns: 2, gutter: 10pt,
  ```python
  from chip_lod.lod import lod_frequentist
  res = lod_frequentist(n=[500, 1000, 2000],
                        k=[3, 1, 3], e=0.001)
  res["called"], res["lod_vaf"]
  ```,
  ```r
  source("r/load.R")
  lod_frequentist(n = c(500, 1000, 2000),
                  k = c(3, 1, 3), e = 0.001)
  # data.frame: n, k, lambda_bg, pvalue,
  #             k_star, called, lod_vaf
  ```,
)

#figure(
  table(
    columns: 7,
    align: right,
    table.header([$n$], [$k$], [$lambda_0$], [p-value], [$k^*$], [called], [LoD VAF]),
    ..f5.map(r => ([#r.n], [#r.k], fmt(r.lambda_bg, d: 3), prob(r.pvalue), [#r.k_star],
      if r.called == "1" [yes] else [no], pct(r.lod_vaf))).flatten(),
  ),
  caption: [Output of `lod_frequentist` for example sites at $e$ = #F.reference_error_rate (`f5_lod_function_example.csv`).],
)

*Methods paragraph* (for a paper; also in `spec/methods_lod_frequentist.txt`):

#block(inset: 10pt, stroke: 0.5pt + luma(160), radius: 4pt, width: 100%, text(size: 9.5pt,
  read("../spec/methods_lod_frequentist.txt")))

// =====================================================================================================
#part[Part II — The Bayesian analysis]

#let pr = B.by_prior
#let pv3 = load("a3_read_table.csv").find(r => int(r.alt_reads) == B.alt_reads).pvalue

= Why another approach?

The p-value answers: _"if there were no variant, how often would we see this many reads?"_ What we usually want
to know is the reverse: _"given that we saw these reads, how likely is it that the variant is real?"_ These are
different questions, and their answers can be very different.

A familiar example is a screening test. A test that rarely gives false alarms can still produce mostly false
alarms when the condition is rare, simply because almost everyone tested does not have it. Variant calling at a
specific hotspot is the same situation. Most people do not carry a given CHIP variant, so most "suspicious" read
counts come from the many people without the variant whose errors happened to pile up.

Bayesian statistics is the bookkeeping that combines two things: *how rare the variant is*, and *how strongly the
reads point to it*. Its result is a probability that the variant is present, which is the number most people
assume the p-value is.

= Bayesian reasoning, step by step

We follow one concrete case throughout: depth *#B.depth*, error rate #B.error_rate, and *#B.alt_reads variant
reads* observed. The frequentist test from Part I #if B.alt_reads >= s.a3.k_star [calls] else [does not call] this
site (p-value #prob(pv3), $k^*$ = #s.a3.k_star). How sure should we be that the variant is real?

== Step 1: what did we expect before looking? (the _prior_)

Imagine #B.cohort people, all sequenced at this position. Before looking at any reads we need a guess for how many
of them carry the variant. This guess is the *prior probability* $pi$. For a specific hotspot in a random person
it is small; for a variant already seen in an earlier sample of the same person it is large. We use two
scenarios throughout:

- *Rare variant*, $pi$ = #pr.at(0).prior: #ppl(pr.at(0).carriers) of the #B.cohort people carry the variant.
- *Agnostic*, $pi$ = #pr.at(1).prior: a coin flip, #ppl(pr.at(1).carriers) carriers.

The prior also says what VAF a variant would have _if_ it is present. We do not know this, so we spread our
belief evenly over orders of magnitude from #pct(s.model.p_min, d: 2) to #pct(s.model.p_max, d: 0): a VAF between
0.1% and 1% is considered as likely as one between 1% and 10%. (Statisticians call this a log-uniform prior.)

== Step 2: what does each explanation predict? (the _likelihood_)

There are two possible explanations for the reads: *errors only*, or *a real variant* (plus errors). Each predicts a
different spread of read counts:

#figure(
  image("../results/figures/fig_b0_likelihood.svg", width: 72%),
  caption: [How often each explanation produces exactly $k$ variant reads at depth #B.depth (`b0_likelihoods.csv`).
    The grey line marks the observed #B.alt_reads reads.],
)

#let b0 = load("b0_likelihoods.csv")
#let k99 = {
  let cum = 0.0
  let out = none
  for r in b0 { cum += float(r.p_error_only); if out == none and cum >= 0.99 { out = r.alt_reads } }
  out
}
- *Errors only* (black) produce at most #k99 reads in over 99% of samples. Exactly #B.alt_reads reads happen in
  #pct(B.m0, d: 2) of variant-free samples.
- *A real variant* (red) produces anything from 0 to many reads, depending on its VAF (dashed lines). Averaged over
  the plausible VAFs from step 1, exactly #B.alt_reads reads happen in #pct(B.m1, d: 2) of carriers.

The ratio of these two numbers is the *Bayes factor*. It summarises the evidence in the reads:

#align(center, keybox[
  Bayes factor = #pct(B.m1, d: 2) / #pct(B.m0, d: 2) = *#fmt(B.bayes_factor, d: 1)*.
  #B.alt_reads reads are #fmt(B.bayes_factor, d: 1) times more likely if the variant is real than if it is not.
])

The Bayes factor depends only on the data and the model, not on $pi$. A value above 1 favours the variant, and the
larger it is, the stronger the evidence. Note that a variant does _not_ always produce many reads: at low VAF it
often produces 0 or 1, which is why the red curve is spread out.

== Step 3: combine the two (the _posterior_)

Now we count, among our #B.cohort imaginary people, *who ends up with exactly #B.alt_reads variant reads*:

#figure(
  table(
    columns: 3,
    align: (left, right, right),
    table.header([], [rare variant ($pi$ = #pr.at(0).prior)], [agnostic ($pi$ = #pr.at(1).prior)]),
    [carriers among #B.cohort], [#ppl(pr.at(0).carriers)], [#ppl(pr.at(1).carriers)],
    [... of whom show #B.alt_reads reads (× #pct(B.m1, d: 2))], [#ppl(pr.at(0).carriers_with_k)], [#ppl(pr.at(1).carriers_with_k)],
    [non-carriers among #B.cohort], [#ppl(pr.at(0).noncarriers)], [#ppl(pr.at(1).noncarriers)],
    [... of whom show #B.alt_reads reads (× #pct(B.m0, d: 2))], [#ppl(pr.at(0).noncarriers_with_k)], [#ppl(pr.at(1).noncarriers_with_k)],
    table.hline(),
    [*share of people with #B.alt_reads reads who are carriers*], [*#pct(pr.at(0).p_h1, d: 1)*], [*#pct(pr.at(1).p_h1, d: 1)*],
  ),
  caption: [The posterior probability as a head count (from `summary.json`, section `bayes_example`).],
)

That last line is the *posterior probability* $P(H_1 | k)$: the probability that the variant is present, given
the reads. ($H_1$ is the hypothesis "variant present"; $H_0$ is "errors only".)

With the rare-variant prior, only #pct(pr.at(0).p_h1, d: 1) of people with #B.alt_reads reads actually carry the
variant. The reads are #fmt(B.bayes_factor, d: 1) times more typical of carriers, but non-carriers outnumber
carriers #calc.round(float(pr.at(0).noncarriers) / float(pr.at(0).carriers)) to 1, so most people with
#B.alt_reads reads are non-carriers with unlucky errors. With the agnostic prior the same reads give
#pct(pr.at(1).p_h1, d: 1). This is the key point of the Bayesian view: *the same reads mean different things
depending on how common the variant is*.

#note[
  *The shortcut with odds.* Bookmakers' odds give the same answer in one line:
  posterior odds = prior odds × Bayes factor. With $pi$ = #pr.at(0).prior the prior odds are
  #fmt(pr.at(0).prior_odds, d: 4) (1 : #calc.round(1 / float(pr.at(0).prior_odds))), so the posterior odds are
  #fmt(pr.at(0).prior_odds, d: 4) × #fmt(B.bayes_factor, d: 1) = #fmt(pr.at(0).posterior_odds, d: 3), i.e. a
  probability of #pct(pr.at(0).p_h1, d: 1). With $pi$ = #pr.at(1).prior: 1 × #fmt(B.bayes_factor, d: 1) gives
  #pct(pr.at(1).p_h1, d: 1).
]

== Step 4: is it CHIP-sized?

A second question uses the same head count: of the carriers with #B.alt_reads reads, how many have a VAF of at
least #pct(s.model.p_thr, d: 0)? With $pi$ = #pr.at(0).prior that is #ppl(pr.at(0).carriers_with_k_vaf_ge_thr) of
#ppl(pr.at(0).carriers_with_k) people. A #pct(s.model.p_thr, d: 0) variant would produce about
#calc.round(B.depth * (s.model.p_thr * (1 - B.error_rate) + (1 - s.model.p_thr) * B.error_rate / 3)) reads at this
depth, not #B.alt_reads. So $P("VAF" >= 2% | k)$ = #prob(pr.at(0).p_vaf_ge_thr): this is almost certainly _not_ a CHIP-sized
clone, even if a variant is present.

== Glossary and decision rules

#table(
  columns: (auto, 1fr),
  [*prior* $pi$], [Probability that the variant is present, before looking at the reads.],
  [*VAF prior*], [What VAF the variant is expected to have if present (here: evenly spread over orders of magnitude between #pct(s.model.p_min, d: 2) and #pct(s.model.p_max, d: 0)).],
  [*likelihood*], [How often an explanation (errors only, or variant) produces the observed read count.],
  [*Bayes factor*], [Likelihood under "variant" divided by likelihood under "errors only": the strength of the evidence in the reads.],
  [*posterior* $P(H_1 | k)$], [Probability that the variant is present, given the reads and the prior.],
  [$P("VAF" >= 2% | k)$], [Probability that the variant is present _and_ at or above the CHIP threshold.],
  [$k_(H 1)$, $k_"thr"$], [Fewest variant reads that push $P(H_1 | k)$, or $P("VAF" >= 2% | k)$, to at least #s.model.tau.],
)

#note[
  *p-value versus posterior.* In our example the p-value is #prob(pv3): significant at $alpha$ = #F.alpha. Yet with a
  rare-variant prior the posterior is only #pct(pr.at(0).p_h1, d: 1). There is no contradiction. The p-value
  says "errors rarely produce this", and that is true. The posterior adds "but variants are rarer still". Report
  both, and state the prior you used.
]

The following three analyses apply this machinery systematically.

= Analysis B1: with few reads, the prior decides

#let a1 = load("a1_prior_sensitivity.csv")
#let ex = s.a1
#let a1-priors = a1.map(r => r.prior).dedup()
#let a1-sub = a1.filter(r => istr(r.depth) == str(ex.depth) and near(r.error_rate, ex.error_rate))
#let a1-ks = a1-sub.map(r => r.alt_reads).dedup()

Take a shallow site: depth #ex.depth with $e$ = #ex.error_rate. The table below shows the p-value and the posterior
$P(H_1 | k)$ for 0 to #a1-ks.last() variant reads under four priors, from very sceptical ($pi$ = #fmt(a1-priors.first(), d: 3))
to agnostic ($pi$ = #fmt(a1-priors.last(), d: 3)).

#figure(
  table(
    columns: 2 + a1-priors.len(),
    align: right,
    table.header(
      table.cell(rowspan: 2)[variant reads $k$], table.cell(rowspan: 2)[p-value],
      table.cell(colspan: a1-priors.len(), align: center)[$P(H_1 | k)$ with prior],
      ..a1-priors.map(p => [$pi$ = #fmt(p, d: 3)])),
    ..a1-ks.map(k => {
      let rs = a1-sub.filter(r => r.alt_reads == k)
      ([#k], prob(rs.first().pvalue), ..a1-priors.map(p => prob(rs.find(r => r.prior == p).p_h1)))
    }).flatten(),
  ),
  caption: [Depth #ex.depth, error rate #ex.error_rate (`a1_prior_sensitivity.csv`).],
)

With #ex.alt_reads variant reads, errors alone are very unlikely (p = #prob(ex.pvalue)). Yet the posterior is
#pct(ex.p_h1_prior_low, d: 2) when variants are rare a priori ($pi$ = #ex.prior_low), and
#pct(ex.p_h1_prior_high, d: 2) under a coin-flip prior. The data are identical; only the prior differs. With a
handful of reads, the evidence (the Bayes factor) is too weak to overrule a strong prior belief.

#figure(
  image("../results/figures/fig_a1_prior.svg", width: 80%),
  caption: [$P(H_1 | k)$ against variant reads at depths 10, 20 and 50 for four priors ($e$ = 0.001).],
)

= Analysis B2: more depth, more certainty

#let a2 = load("a2_depth_curve_analytic.csv")
#let a2s = load("a2_depth_curve_sim.csv")
#let a2-priors = a2.map(r => r.prior).dedup()
#let a2-depths = a2.map(r => r.depth).dedup()

Now fix a true variant at VAF #pct(s.a2.true_vaf, d: 0), the CHIP threshold, with $e$ = #s.a2.error_rate, and increase
the depth. A real sample gives one particular read count; another sample of the same variant might give a few more
or fewer. So for each depth we report the *expected* posterior: the average of $P(H_1 | k)$ over all the read counts
such a variant could produce, weighted by how often each occurs.

#figure(
  table(
    columns: 3 + 2 * a2-priors.len(),
    align: right,
    table.header(
      table.cell(rowspan: 2)[depth], table.cell(rowspan: 2)[$k^*$], table.cell(rowspan: 2)[LoD VAF],
      table.cell(colspan: a2-priors.len(), align: center)[expected $P(H_1 | k)$],
      table.cell(colspan: a2-priors.len(), align: center)[expected $P("VAF" >= #pct(s.model.p_thr, d: 0) | k)$],
      ..(a2-priors + a2-priors).map(p => [$pi$ = #fmt(p)])),
    ..a2-depths.map(d => {
      let rs = a2.filter(r => r.depth == d)
      ([#d], [#rs.first().k_star], pct(rs.first().lod_vaf),
        ..a2-priors.map(p => prob(rs.find(r => r.prior == p).expected_p_h1)),
        ..a2-priors.map(p => prob(rs.find(r => r.prior == p).expected_p_vaf_ge_thr)))
    }).flatten(),
  ),
  caption: [Expected posteriors for a true VAF of #pct(s.a2.true_vaf, d: 0); $k^*$ and LoD at $alpha$ = #s.model.alpha
    for comparison (`a2_depth_curve_analytic.csv`).],
)

*Reading the table.* With more reads, the evidence grows until it overwhelms any reasonable prior. The expected
$P(H_1)$ first reaches #s.model.tau at #dfmt(s.a2.by_prior.at(0).depth_expected_p_h1_ge_tau) with $pi$ =
#s.a2.by_prior.at(0).prior, and #dfmt(s.a2.by_prior.at(1).depth_expected_p_h1_ge_tau) with $pi$ =
#s.a2.by_prior.at(1).prior: a sceptical prior needs more evidence, so more depth.

The last columns show a subtle point. Because the true VAF here is _exactly_ the 2% threshold, the probability
that the VAF is "at least 2%" does not go to 1. It levels off near 0.5 (#prob(a2.last().expected_p_vaf_ge_thr) at
#a2.last().depth reads), because a variant sitting on the line is equally likely to look just above or just below
it. Concluding that a variant is _over_ the CHIP threshold needs a true VAF clearly above it.

#figure(
  image("../results/figures/fig_a2_depth.svg", width: 72%),
  caption: [Expected $P(H_1)$ against depth (lines) with simulated sites: median (points) and 10th–90th percentile
    band at VAF 2%; dashed lines show the median at VAF 0 (no variant). Dotted line: 0.95.],
)

*Simulated check.* We simulated #istr(a2s.first().n_sites) sites per depth at VAF 2% and at VAF 0
(`a2_depth_curve_sim.csv`). The analytic curve is an _average_, so it is compared with the whole simulated spread
(shaded band). The band closes in on 1 at the same depths where the curve does. Among the #s.a2.n_null_sites
sites with _no_ variant, the frequentist test called #pct(s.a2.sim_fp_rate_freq) as positive, while the rule
"$P(H_1 | k) >= #s.model.tau$" called #pct(s.a2.by_prior.at(0).sim_fp_rate_bayes) ($pi$ = #s.a2.by_prior.at(0).prior)
and #pct(s.a2.by_prior.at(1).sim_fp_rate_bayes) ($pi$ = #s.a2.by_prior.at(1).prior).

= Analysis B3: at depth #s.a3.depth, how many reads are enough?

#let a3r = load("a3_read_table.csv")
#let a3t = load("a3_thresholds.csv")
#let a3s = load("a3_sim_sensitivity.csv")
#let sim3-depths = csv(res + "sim_a3.csv").slice(1).map(r => int(r.at(2)))
#let deep = a3t.filter(r => near(r.alpha, s.a3.alpha)).last()
#let a3-priors = a3r.map(r => r.prior).dedup()
#let a3-ks = a3r.map(r => r.alt_reads).dedup()
#let show-k = a3-ks.filter(k => float(k) <= 12 or calc.rem(calc.round(float(k)), 5) == 0)

This is the everyday case: a panel with about #s.a3.depth reads at the hotspot and $e$ = #s.a3.error_rate.

#figure(
  table(
    columns: 2 + 2 * a3-priors.len(),
    align: right,
    table.header(
      table.cell(rowspan: 2)[variant reads $k$], table.cell(rowspan: 2)[p-value],
      table.cell(colspan: a3-priors.len(), align: center)[$P(H_1 | k)$],
      table.cell(colspan: a3-priors.len(), align: center)[$P("VAF" >= #pct(s.model.p_thr, d: 0) | k)$],
      ..(a3-priors + a3-priors).map(p => [$pi$ = #fmt(p)])),
    ..show-k.map(k => {
      let rs = a3r.filter(r => r.alt_reads == k)
      ([#k], prob(rs.first().pvalue),
        ..a3-priors.map(p => prob(rs.find(r => r.prior == p).p_h1)),
        ..a3-priors.map(p => prob(rs.find(r => r.prior == p).p_vaf_ge_thr)))
    }).flatten(),
  ),
  caption: [Evidence for each variant read count at depth #s.a3.depth (selected rows of `a3_read_table.csv`, which has
    every $k$ from 0 to #a3-ks.last()).],
)

#figure(
  image("../results/figures/fig_a3_reads.svg", width: 70%),
  caption: [Posteriors against variant reads at depth #s.a3.depth. Solid: $P(H_1 | k)$; dashed: $P("VAF" >= 2% | k)$
    (the two priors overlap). Dotted vertical line: the frequentist $k^*$.],
)

*Three thresholds, three questions.*
- *Is it more than error?* (frequentist) At least #s.a3.k_star variant reads ($k^*$, $alpha$ = #s.a3.alpha). Variants
  at #pct(s.a3.lod_vaf) VAF or higher clear this bar in #pct(s.model.power, d: 0) of samples.
- *Do I believe it is there?* (Bayesian) At least #s.a3.by_prior.at(0).k_h1 reads if variants are rare
  ($pi$ = #s.a3.by_prior.at(0).prior), or #s.a3.by_prior.at(1).k_h1 with a coin-flip prior ($k_(H 1)$).
- *Is it CHIP-sized?* (Bayesian) At least #s.a3.by_prior.at(0).k_thr reads for 95% confidence that the VAF is at least
  #pct(s.model.p_thr, d: 0) ($k_"thr"$) with $pi$ = #s.a3.by_prior.at(0).prior, and #s.a3.by_prior.at(1).k_thr with
  $pi$ = #s.a3.by_prior.at(1).prior#if s.a3.by_prior.at(0).k_thr == s.a3.by_prior.at(1).k_thr [: here the prior no
  longer matters, because with this many reads the evidence dominates] else [].

#figure(
  table(
    columns: 7,
    align: right,
    table.header([depth], [$alpha$], [prior $pi$], [$k^*$], [LoD VAF], [$k_(H 1)$], [$k_"thr"$]),
    ..a3t.map(r => ([#r.depth], [#fmt(r.alpha, d: 3)], [#fmt(r.prior)], [#r.k_star], pct(r.lod_vaf),
      [#r.k_h1], [#r.k_thr])).flatten(),
  ),
  caption: [Thresholds for depths around #s.a3.depth (`a3_thresholds.csv`). $k^*$ and the LoD do not depend on the prior;
    $k_(H 1)$ and $k_"thr"$ do not depend on $alpha$.],
)

*Simulated sensitivity.* We simulated sites whose depth varies uniformly between #calc.min(..sim3-depths) and
#calc.max(..sim3-depths) at several true VAFs. Each site was called with the frequentist rule (its own $k^*$) and with
the Bayesian rule ($P(H_1 | k) >= #s.model.tau$, $pi$ = #s.a3.sim_prior).

#figure(
  table(
    columns: 6,
    align: right,
    table.header([true VAF], [sites], [mean depth], [frequentist call rate], [Bayesian call rate],
      [analytic power at #s.a3.depth]),
    ..a3s.map(r => (pct(r.true_vaf), [#r.n_sites], fmt(r.mean_depth, d: 0), pct(r.call_rate_freq, d: 1),
      pct(r.call_rate_bayes, d: 1), pct(r.analytic_power_at_1000, d: 1))).flatten(),
  ),
  caption: [Call rates in simulated data (`a3_sim_sensitivity.csv`). The VAF 0 row is the false-positive rate.],
)

#figure(
  image("../results/figures/fig_a3_sensitivity.svg", width: 70%),
  caption: [Call rate against true VAF in simulated sites (99% Wilson intervals), with the analytic power curve.],
)

With no variant, the frequentist rule called #pct(s.a3.sim_fp_rate_freq, d: 1) of sites and the Bayesian rule
#pct(s.a3.sim_fp_rate_bayes, d: 1). The Bayesian rule with a sceptical prior gives up sensitivity at the lowest
VAFs in exchange for far fewer false calls. The analytic power assumes exactly #s.a3.depth reads; deeper simulated
sites need more reads ($k^*$ = #deep.k_star at depth #deep.depth), which is why the frequentist call rates sit
slightly below it.

= Which approach should I use?

#table(
  columns: (auto, 1fr, 1fr),
  table.header([], [*frequentist (Part I)*], [*Bayesian (Part II)*]),
  [question answered], [Are these reads more than errors would produce?], [How likely is the variant to be real (and CHIP-sized)?],
  [needs], [error rate, $alpha$], [error rate, prior $pi$, VAF prior],
  [main output], [$k^*$, p-value, limit of detection], [posterior probability, $k_(H 1)$, $k_"thr"$],
  [strength], [simple, assay-level, standard in method papers], [directly interpretable, uses what is known about the variant],
  [weakness], [ignores how rare the variant is], [depends on the prior, which must be justified],
)

For assay validation and method sections, the frequentist LoD (Part I, with `lod_frequentist`) is the standard
measure. For interpreting an individual sample, the posterior adds what the p-value cannot: how much to believe
the call given how common the variant is.

= Reproducibility

The package has two independent implementations, Python and R, written to produce *byte-identical* result files.
Both use the same hand-written random number generator, explicit left-to-right sums and the same number
formatting. `scripts/compare_outputs.sh` compares the SHA-256 hash of every file.

#let hfiles = hc.files.keys()
#let n-match = hfiles.filter(f => hc.files.at(f).match).len()
#block(fill: if hc.all_match { rgb("#e6f2e6") } else { rgb("#f8e0e0") }, inset: 8pt, radius: 4pt, width: 100%)[
  *#if hc.all_match [Python and R agree:] else [Python and R DISAGREE:]* #n-match of #hfiles.len() result files
  have identical SHA-256 hashes.
]

#table(
  columns: 3,
  table.header([file], [SHA-256 (first 16 hex digits)], [match]),
  ..hfiles.map(f => {
    let e = hc.files.at(f)
    (raw(f), raw(if e.python == none { "missing" } else { e.python.slice(0, 16) }), if e.match [yes] else [*no*])
  }).flatten(),
)

Every number in this report is read from `results/python/` when the document is compiled, so rerunning `make all`
updates everything. The model is specified in `spec/MODEL.md`, and every judgement call is listed in
`spec/DECISIONS.md`.

= Limitations

- *Poisson approximation.* Read counts are really binomial, and the Poisson is slightly too wide. The simulations
  show the difference is negligible at these depths and VAFs.
- *No overdispersion.* Real error counts vary more from sample to sample than a Poisson allows (batch, library and
  position effects). A beta-binomial or negative-binomial background would give larger, more honest $k^*$ values.
- *The error rate is site- and context-specific.* A single $e$ = #F.reference_error_rate is a placeholder. Error rates
  depend on the base change, the sequence context and the assay; Part I shows how strongly this matters. In practice,
  estimate the background at each hotspot from a *panel of normals*.
- *PCR duplicates and UMIs.* The model treats reads as independent. Duplicates from the same molecule are not
  independent and make a variant look better supported than it is. Count unique molecules (UMI families), not raw
  reads; UMI consensus also lowers $e$ considerably.
- *Diploid, heterozygous assumption.* VAF = cell fraction / 2 assumes a heterozygous variant on a normal copy number.
  Loss of heterozygosity, copy-number changes or a homozygous variant break this link.
- *The prior is a choice.* The posterior is only as good as $pi$ and the VAF prior. Analysis B1 shows that with few
  reads it can dominate. Use a prior grounded in population frequencies for the specific hotspot, and report it.
