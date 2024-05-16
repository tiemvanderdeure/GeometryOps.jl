import{_ as e,c as i,o as l,a6 as t}from"./chunks/framework.BZUDmi18.js";const m=JSON.parse('{"title":"","description":"","frontmatter":{},"headers":[],"relativePath":"call_notes.md","filePath":"call_notes.md","lastUpdated":null}'),o={name:"call_notes.md"},a=t('<h2 id="20th-April,-2024" tabindex="-1">20th April, 2024 <a class="header-anchor" href="#20th-April,-2024" aria-label="Permalink to &quot;20th April, 2024 {#20th-April,-2024}&quot;">​</a></h2><p>See <a href="https://github.com/JuliaGeo/GeometryOps.jl/issues/114" target="_blank" rel="noreferrer">GeometryOps#114</a>.</p><ul><li><p>[ ] Exact predicates can be defined for lower-level, more atomic predicates within GeometryOps.</p></li><li><p>[ ] Add Shewchuck&#39;s adaptive math as a stage for exact predicates.</p><ul><li>[ ] We can add this by simply adding a function which creates and does the error check for MultiFloat numbers, in the vein of <a href="https://github.com/lairez/ExactPredicates.jl/blob/a9dce0334ee62104a14930fe206764bf802c23db/src/Codegen.jl#L468C1-L471C13" target="_blank" rel="noreferrer">the code to generate the main function of a predicate</a> from ExactPredicates.</li></ul></li><li><p>[x] @skygering to write docstrings for the predicates</p></li></ul><h2 id="29th-Feb,-2024" tabindex="-1">29th Feb, 2024 <a class="header-anchor" href="#29th-Feb,-2024" aria-label="Permalink to &quot;29th Feb, 2024 {#29th-Feb,-2024}&quot;">​</a></h2><h3 id="To-do" tabindex="-1">To do <a class="header-anchor" href="#To-do" aria-label="Permalink to &quot;To do {#To-do}&quot;">​</a></h3><ul><li><p>[ ] Finish clipping degeneracies</p></li><li><p>[ ] Fix cross &amp; overlap functions</p></li><li><p>[x] Benchmarks to show why things you couldn&#39;t concieve of in R are doable in Julia</p></li><li><p>[x] profile functions for exponential improvements</p></li><li><p>[ ] A list of projects people can work on...the beauty here is that each function is kind of self-contained so it&#39;s an undergrad level project</p></li><li><p>[ ] Doc improvements</p><ul><li><p>more</p></li><li><p>benchmarks page</p></li></ul></li><li><p>Methods to validate and fix geometry</p><ul><li><p>[ ] Polygons and LinearRings:</p><ul><li><p>[ ] self-intersection</p></li><li><p>[ ] holes are actually within the polygon</p></li><li><p>[ ] Polygon exteriors must be counterclockwise, holes clockwise.</p></li><li><p>[ ] length of all rings &gt; 4</p></li><li><p>[ ] repeated last point</p></li></ul></li><li><p>[ ] LineStrings: NaN/Inf points</p></li><li><p>[x] Fix linear rings at some point to make sure the ring is closed, i.e., <code>points[end] == points[begin]</code></p></li><li></li><li></li></ul></li><li><p>Tests</p><ul><li><p>[x] Simplify functions</p></li><li><p>[x] Polygonize</p></li><li><p>Barycentric tests for <code>n_vertices &gt; 4</code></p></li></ul></li></ul><h3 id="Done" tabindex="-1">Done <a class="header-anchor" href="#Done" aria-label="Permalink to &quot;Done {#Done}&quot;">​</a></h3><ul><li><p>Rename <code>bools.jl</code> to something more relevant to the actual code -&gt; <code>orientation.jl</code></p></li><li><p>Doc improvements:</p><ul><li>organise sections</li></ul></li></ul>',8),n=[a];function r(s,p,c,d,h,u){return l(),i("div",null,n)}const g=e(o,[["render",r]]);export{m as __pageData,g as default};
