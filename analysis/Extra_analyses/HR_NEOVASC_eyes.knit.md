---
title: "Hyper-Reflectivity and Neovascularization in LMRI DB"
author: "Christian Anderson"
date: "2025-02-25"
output:
  html_document:
    df_print: paged
  always_allow_html: true
  pdf_document:
    keep_tex: true
---


















### Q: What percent of MacTel patients have outer retina hyper reflectivity (ORHR) and choroidal neovascularization (CNV)?
ORHR and CNV are not variables in the database, but HR_TEMPORAL (Hyper reflectivity on the temporal side) and NEOVASC (Neovascularization) are. 

__4.9__ % or eyes have both HR_TEMPORAL and NEOVASC while __31.7%__ only have HR_TEMPORAL and __1.2%__ only NEOVASC. These percentages are found in Table 2.

If you are interested in seeing a breakdown of phenotype presence on the level of patients and eyes, please reference the tables below. 

#### Tables
 
 
 

```{=html}
<div id="quljooctka" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#quljooctka table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#quljooctka thead, #quljooctka tbody, #quljooctka tfoot, #quljooctka tr, #quljooctka td, #quljooctka th {
  border-style: none;
}

#quljooctka p {
  margin: 0;
  padding: 0;
}

#quljooctka .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 600px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#quljooctka .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#quljooctka .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#quljooctka .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#quljooctka .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#quljooctka .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#quljooctka .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#quljooctka .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#quljooctka .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#quljooctka .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#quljooctka .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#quljooctka .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#quljooctka .gt_spanner_row {
  border-bottom-style: hidden;
}

#quljooctka .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#quljooctka .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#quljooctka .gt_from_md > :first-child {
  margin-top: 0;
}

#quljooctka .gt_from_md > :last-child {
  margin-bottom: 0;
}

#quljooctka .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#quljooctka .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#quljooctka .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#quljooctka .gt_row_group_first td {
  border-top-width: 2px;
}

#quljooctka .gt_row_group_first th {
  border-top-width: 2px;
}

#quljooctka .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#quljooctka .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#quljooctka .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#quljooctka .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#quljooctka .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#quljooctka .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#quljooctka .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#quljooctka .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#quljooctka .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#quljooctka .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#quljooctka .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#quljooctka .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#quljooctka .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#quljooctka .gt_left {
  text-align: left;
}

#quljooctka .gt_center {
  text-align: center;
}

#quljooctka .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#quljooctka .gt_font_normal {
  font-weight: normal;
}

#quljooctka .gt_font_bold {
  font-weight: bold;
}

#quljooctka .gt_font_italic {
  font-style: italic;
}

#quljooctka .gt_super {
  font-size: 65%;
}

#quljooctka .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#quljooctka .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#quljooctka .gt_indent_1 {
  text-indent: 5px;
}

#quljooctka .gt_indent_2 {
  text-indent: 10px;
}

#quljooctka .gt_indent_3 {
  text-indent: 15px;
}

#quljooctka .gt_indent_4 {
  text-indent: 20px;
}

#quljooctka .gt_indent_5 {
  text-indent: 25px;
}

#quljooctka .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#quljooctka div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="5" class="gt_heading gt_title gt_font_normal" style>Eyes that develop NEOVASC and/or HR_TEMPORAL</td>
    </tr>
    <tr class="gt_heading">
      <td colspan="5" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>0: not Present; 1: Questionable; 2: Present</td>
    </tr>
    <tr class="gt_col_headings gt_spanner_row">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="2" scope="colgroup" id="a::stub"></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="3" scope="colgroup" id="HR_TEMPORAL">
        <div class="gt_column_spanner">HR_TEMPORAL</div>
      </th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="a0">0</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="a1">1</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="a2">2</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_row_group_first"><td headers="NEOVASC stub_2_1 stub_1" rowspan="3" class="gt_row gt_left gt_stub_row_group">NEOVASC</td>
<th id="stub_2_1" scope="row" class="gt_row gt_right gt_stub">0</th>
<td headers="NEOVASC stub_2_1 0" class="gt_row gt_right">4535</td>
<td headers="NEOVASC stub_2_1 1" class="gt_row gt_right">1039</td>
<td headers="NEOVASC stub_2_1 2" class="gt_row gt_right">1266</td></tr>
    <tr><th id="0_2" scope="row" class="gt_row gt_right gt_stub">1</th>
<td headers="NEOVASC 0_2 0" class="gt_row gt_right">63</td>
<td headers="NEOVASC 0_2 1" class="gt_row gt_right">47</td>
<td headers="NEOVASC 0_2 2" class="gt_row gt_right">120</td></tr>
    <tr><th id="0_3" scope="row" class="gt_row gt_right gt_stub">2</th>
<td headers="NEOVASC 0_3 0" class="gt_row gt_right">21</td>
<td headers="NEOVASC 0_3 1" class="gt_row gt_right">20</td>
<td headers="NEOVASC 0_3 2" class="gt_row gt_right">174</td></tr>
  </tbody>
  <tfoot class="gt_sourcenotes">
    <tr>
      <td class="gt_sourcenote" colspan="5">X^2 test: p ~ 1.061e-138; N = 7285</td>
    </tr>
    <tr>
      <td class="gt_sourcenote" colspan="5">Note: These are counts of unique eyes that ever develop phenotype, not counts of eye-visits</td>
    </tr>
  </tfoot>
  
</table>
</div>
```

__Table 1__: Counts of eyes for every combination of HR_TEMPORAL and NEOVASC. A chi-square test was performed on these values to test whether the presence of phenotype X depends on phenotype Y. 
 
-----------------------------------
 

```{=html}
<div id="yfuintgixt" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#yfuintgixt table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#yfuintgixt thead, #yfuintgixt tbody, #yfuintgixt tfoot, #yfuintgixt tr, #yfuintgixt td, #yfuintgixt th {
  border-style: none;
}

#yfuintgixt p {
  margin: 0;
  padding: 0;
}

#yfuintgixt .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 600px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#yfuintgixt .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#yfuintgixt .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#yfuintgixt .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#yfuintgixt .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#yfuintgixt .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#yfuintgixt .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#yfuintgixt .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#yfuintgixt .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#yfuintgixt .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#yfuintgixt .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#yfuintgixt .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#yfuintgixt .gt_spanner_row {
  border-bottom-style: hidden;
}

#yfuintgixt .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#yfuintgixt .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#yfuintgixt .gt_from_md > :first-child {
  margin-top: 0;
}

#yfuintgixt .gt_from_md > :last-child {
  margin-bottom: 0;
}

#yfuintgixt .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#yfuintgixt .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#yfuintgixt .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#yfuintgixt .gt_row_group_first td {
  border-top-width: 2px;
}

#yfuintgixt .gt_row_group_first th {
  border-top-width: 2px;
}

#yfuintgixt .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#yfuintgixt .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#yfuintgixt .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#yfuintgixt .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#yfuintgixt .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#yfuintgixt .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#yfuintgixt .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#yfuintgixt .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#yfuintgixt .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#yfuintgixt .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#yfuintgixt .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#yfuintgixt .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#yfuintgixt .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#yfuintgixt .gt_left {
  text-align: left;
}

#yfuintgixt .gt_center {
  text-align: center;
}

#yfuintgixt .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#yfuintgixt .gt_font_normal {
  font-weight: normal;
}

#yfuintgixt .gt_font_bold {
  font-weight: bold;
}

#yfuintgixt .gt_font_italic {
  font-style: italic;
}

#yfuintgixt .gt_super {
  font-size: 65%;
}

#yfuintgixt .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#yfuintgixt .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#yfuintgixt .gt_indent_1 {
  text-indent: 5px;
}

#yfuintgixt .gt_indent_2 {
  text-indent: 10px;
}

#yfuintgixt .gt_indent_3 {
  text-indent: 15px;
}

#yfuintgixt .gt_indent_4 {
  text-indent: 20px;
}

#yfuintgixt .gt_indent_5 {
  text-indent: 25px;
}

#yfuintgixt .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#yfuintgixt div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="5" class="gt_heading gt_title gt_font_normal" style>Ratio of eyes that develop NEOVASC and/or HR_TEMPORAL</td>
    </tr>
    <tr class="gt_heading">
      <td colspan="5" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>0: not Present; 1: Questionable; 2: Present</td>
    </tr>
    <tr class="gt_col_headings gt_spanner_row">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="2" scope="colgroup" id="a::stub"></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="3" scope="colgroup" id="HR_TEMPORAL">
        <div class="gt_column_spanner">HR_TEMPORAL</div>
      </th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="a0">0</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="a1">1</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="a2">2</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_row_group_first"><td headers="NEOVASC stub_2_1 stub_1" rowspan="3" class="gt_row gt_left gt_stub_row_group">NEOVASC</td>
<th id="stub_2_1" scope="row" class="gt_row gt_right gt_stub">0</th>
<td headers="NEOVASC stub_2_1 0" class="gt_row gt_right">0.623</td>
<td headers="NEOVASC stub_2_1 1" class="gt_row gt_right">0.143</td>
<td headers="NEOVASC stub_2_1 2" class="gt_row gt_right">0.174</td></tr>
    <tr><th id="0_2" scope="row" class="gt_row gt_right gt_stub">1</th>
<td headers="NEOVASC 0_2 0" class="gt_row gt_right">0.009</td>
<td headers="NEOVASC 0_2 1" class="gt_row gt_right">0.006</td>
<td headers="NEOVASC 0_2 2" class="gt_row gt_right">0.016</td></tr>
    <tr><th id="0_3" scope="row" class="gt_row gt_right gt_stub">2</th>
<td headers="NEOVASC 0_3 0" class="gt_row gt_right">0.003</td>
<td headers="NEOVASC 0_3 1" class="gt_row gt_right">0.003</td>
<td headers="NEOVASC 0_3 2" class="gt_row gt_right">0.024</td></tr>
  </tbody>
  <tfoot class="gt_sourcenotes">
    <tr>
      <td class="gt_sourcenote" colspan="5">Ratio develop both: 0.049000 (1,1,+1,2+2,2+2,1)</td>
    </tr>
    <tr>
      <td class="gt_sourcenote" colspan="5">Ratio develop HR_TEMPORAL only: 0.317000 (0,1+0,2)</td>
    </tr>
    <tr>
      <td class="gt_sourcenote" colspan="5">Ratio develop NEOVASC only: 0.012000 (1,0+2,0)</td>
    </tr>
  </tfoot>
  
</table>
</div>
```

__Table 2__: Every combination of HR_TEMPORAL and NEOVASC as a ratio of total eyes ('N' in Table 1). The percentages of eyes that develop both and each individually were calculated by summing their respective cells in the table. 
 
-----------------------------------
 

```{=html}
<div id="fpxkhtmfcu" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#fpxkhtmfcu table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#fpxkhtmfcu thead, #fpxkhtmfcu tbody, #fpxkhtmfcu tfoot, #fpxkhtmfcu tr, #fpxkhtmfcu td, #fpxkhtmfcu th {
  border-style: none;
}

#fpxkhtmfcu p {
  margin: 0;
  padding: 0;
}

#fpxkhtmfcu .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 600px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#fpxkhtmfcu .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#fpxkhtmfcu .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#fpxkhtmfcu .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#fpxkhtmfcu .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#fpxkhtmfcu .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#fpxkhtmfcu .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#fpxkhtmfcu .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#fpxkhtmfcu .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#fpxkhtmfcu .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#fpxkhtmfcu .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#fpxkhtmfcu .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#fpxkhtmfcu .gt_spanner_row {
  border-bottom-style: hidden;
}

#fpxkhtmfcu .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#fpxkhtmfcu .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#fpxkhtmfcu .gt_from_md > :first-child {
  margin-top: 0;
}

#fpxkhtmfcu .gt_from_md > :last-child {
  margin-bottom: 0;
}

#fpxkhtmfcu .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#fpxkhtmfcu .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#fpxkhtmfcu .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#fpxkhtmfcu .gt_row_group_first td {
  border-top-width: 2px;
}

#fpxkhtmfcu .gt_row_group_first th {
  border-top-width: 2px;
}

#fpxkhtmfcu .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#fpxkhtmfcu .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#fpxkhtmfcu .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#fpxkhtmfcu .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#fpxkhtmfcu .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#fpxkhtmfcu .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#fpxkhtmfcu .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#fpxkhtmfcu .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#fpxkhtmfcu .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#fpxkhtmfcu .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#fpxkhtmfcu .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#fpxkhtmfcu .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#fpxkhtmfcu .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#fpxkhtmfcu .gt_left {
  text-align: left;
}

#fpxkhtmfcu .gt_center {
  text-align: center;
}

#fpxkhtmfcu .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#fpxkhtmfcu .gt_font_normal {
  font-weight: normal;
}

#fpxkhtmfcu .gt_font_bold {
  font-weight: bold;
}

#fpxkhtmfcu .gt_font_italic {
  font-style: italic;
}

#fpxkhtmfcu .gt_super {
  font-size: 65%;
}

#fpxkhtmfcu .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#fpxkhtmfcu .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#fpxkhtmfcu .gt_indent_1 {
  text-indent: 5px;
}

#fpxkhtmfcu .gt_indent_2 {
  text-indent: 10px;
}

#fpxkhtmfcu .gt_indent_3 {
  text-indent: 15px;
}

#fpxkhtmfcu .gt_indent_4 {
  text-indent: 20px;
}

#fpxkhtmfcu .gt_indent_5 {
  text-indent: 25px;
}

#fpxkhtmfcu .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#fpxkhtmfcu div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="5" class="gt_heading gt_title gt_font_normal" style>Patients that develop NEOVASC and/or HR_TEMPORAL</td>
    </tr>
    <tr class="gt_heading">
      <td colspan="5" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>0: not Present; 1: Questionable; 2: Present</td>
    </tr>
    <tr class="gt_col_headings gt_spanner_row">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="2" scope="colgroup" id="a::stub"></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="3" scope="colgroup" id="HR_TEMPORAL">
        <div class="gt_column_spanner">HR_TEMPORAL</div>
      </th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="a0">0</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="a1">1</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="a2">2</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_row_group_first"><td headers="NEOVASC stub_2_1 stub_1" rowspan="3" class="gt_row gt_left gt_stub_row_group">NEOVASC</td>
<th id="stub_2_1" scope="row" class="gt_row gt_right gt_stub">0</th>
<td headers="NEOVASC stub_2_1 0" class="gt_row gt_right">1869</td>
<td headers="NEOVASC stub_2_1 1" class="gt_row gt_right">614</td>
<td headers="NEOVASC stub_2_1 2" class="gt_row gt_right">810</td></tr>
    <tr><th id="0_2" scope="row" class="gt_row gt_right gt_stub">1</th>
<td headers="NEOVASC 0_2 0" class="gt_row gt_right">36</td>
<td headers="NEOVASC 0_2 1" class="gt_row gt_right">40</td>
<td headers="NEOVASC 0_2 2" class="gt_row gt_right">106</td></tr>
    <tr><th id="0_3" scope="row" class="gt_row gt_right gt_stub">2</th>
<td headers="NEOVASC 0_3 0" class="gt_row gt_right">14</td>
<td headers="NEOVASC 0_3 1" class="gt_row gt_right">14</td>
<td headers="NEOVASC 0_3 2" class="gt_row gt_right">158</td></tr>
  </tbody>
  <tfoot class="gt_sourcenotes">
    <tr>
      <td class="gt_sourcenote" colspan="5">X^2 test: p ~ 9.332e-87; N = 3661</td>
    </tr>
    <tr>
      <td class="gt_sourcenote" colspan="5">Note: These are counts of unique patients that ever develop phenotype, not counts of patient-visits</td>
    </tr>
  </tfoot>
  
</table>
</div>
```

__Table 3__: Counts of patients for every combination of HR_TEMPORAL and NEOVASC. A chi-square test was performed on these values to test whether the presence of phenotype X depends on phenotype Y.
 
-----------------------------------
 

```{=html}
<div id="rqjmbsehbm" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#rqjmbsehbm table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#rqjmbsehbm thead, #rqjmbsehbm tbody, #rqjmbsehbm tfoot, #rqjmbsehbm tr, #rqjmbsehbm td, #rqjmbsehbm th {
  border-style: none;
}

#rqjmbsehbm p {
  margin: 0;
  padding: 0;
}

#rqjmbsehbm .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 600px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#rqjmbsehbm .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#rqjmbsehbm .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#rqjmbsehbm .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#rqjmbsehbm .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#rqjmbsehbm .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rqjmbsehbm .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#rqjmbsehbm .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#rqjmbsehbm .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#rqjmbsehbm .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#rqjmbsehbm .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#rqjmbsehbm .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#rqjmbsehbm .gt_spanner_row {
  border-bottom-style: hidden;
}

#rqjmbsehbm .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#rqjmbsehbm .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#rqjmbsehbm .gt_from_md > :first-child {
  margin-top: 0;
}

#rqjmbsehbm .gt_from_md > :last-child {
  margin-bottom: 0;
}

#rqjmbsehbm .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#rqjmbsehbm .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#rqjmbsehbm .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#rqjmbsehbm .gt_row_group_first td {
  border-top-width: 2px;
}

#rqjmbsehbm .gt_row_group_first th {
  border-top-width: 2px;
}

#rqjmbsehbm .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#rqjmbsehbm .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#rqjmbsehbm .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#rqjmbsehbm .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rqjmbsehbm .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#rqjmbsehbm .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#rqjmbsehbm .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#rqjmbsehbm .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#rqjmbsehbm .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rqjmbsehbm .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#rqjmbsehbm .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#rqjmbsehbm .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#rqjmbsehbm .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#rqjmbsehbm .gt_left {
  text-align: left;
}

#rqjmbsehbm .gt_center {
  text-align: center;
}

#rqjmbsehbm .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#rqjmbsehbm .gt_font_normal {
  font-weight: normal;
}

#rqjmbsehbm .gt_font_bold {
  font-weight: bold;
}

#rqjmbsehbm .gt_font_italic {
  font-style: italic;
}

#rqjmbsehbm .gt_super {
  font-size: 65%;
}

#rqjmbsehbm .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#rqjmbsehbm .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#rqjmbsehbm .gt_indent_1 {
  text-indent: 5px;
}

#rqjmbsehbm .gt_indent_2 {
  text-indent: 10px;
}

#rqjmbsehbm .gt_indent_3 {
  text-indent: 15px;
}

#rqjmbsehbm .gt_indent_4 {
  text-indent: 20px;
}

#rqjmbsehbm .gt_indent_5 {
  text-indent: 25px;
}

#rqjmbsehbm .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#rqjmbsehbm div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="5" class="gt_heading gt_title gt_font_normal" style>Ratio of patients that develop NEOVASC and/or HR_TEMPORAL</td>
    </tr>
    <tr class="gt_heading">
      <td colspan="5" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>0: not Present; 1: Questionable; 2: Present</td>
    </tr>
    <tr class="gt_col_headings gt_spanner_row">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="2" scope="colgroup" id="a::stub"></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="3" scope="colgroup" id="HR_TEMPORAL">
        <div class="gt_column_spanner">HR_TEMPORAL</div>
      </th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="a0">0</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="a1">1</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="a2">2</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_row_group_first"><td headers="NEOVASC stub_2_1 stub_1" rowspan="3" class="gt_row gt_left gt_stub_row_group">NEOVASC</td>
<th id="stub_2_1" scope="row" class="gt_row gt_right gt_stub">0</th>
<td headers="NEOVASC stub_2_1 0" class="gt_row gt_right">0.511</td>
<td headers="NEOVASC stub_2_1 1" class="gt_row gt_right">0.168</td>
<td headers="NEOVASC stub_2_1 2" class="gt_row gt_right">0.221</td></tr>
    <tr><th id="0_2" scope="row" class="gt_row gt_right gt_stub">1</th>
<td headers="NEOVASC 0_2 0" class="gt_row gt_right">0.010</td>
<td headers="NEOVASC 0_2 1" class="gt_row gt_right">0.011</td>
<td headers="NEOVASC 0_2 2" class="gt_row gt_right">0.029</td></tr>
    <tr><th id="0_3" scope="row" class="gt_row gt_right gt_stub">2</th>
<td headers="NEOVASC 0_3 0" class="gt_row gt_right">0.004</td>
<td headers="NEOVASC 0_3 1" class="gt_row gt_right">0.004</td>
<td headers="NEOVASC 0_3 2" class="gt_row gt_right">0.043</td></tr>
  </tbody>
  <tfoot class="gt_sourcenotes">
    <tr>
      <td class="gt_sourcenote" colspan="5">Ratio develop both: 0.087000 (1,1,+1,2+2,2+2,1)</td>
    </tr>
    <tr>
      <td class="gt_sourcenote" colspan="5">Ratio develop HR_TEMPORAL only: 0.389000 (0,1+0,2)</td>
    </tr>
    <tr>
      <td class="gt_sourcenote" colspan="5">Ratio develop NEOVASC only: 0.014000 (1,0+2,0)</td>
    </tr>
  </tfoot>
  
</table>
</div>
```

__Table 4__: Every combination of HR_TEMPORAL and NEOVASC as a ratio of total patients ('N' in Table 3). The percentages of patients that develop both and each individually were calculated by summing their respective cells in the table. 

