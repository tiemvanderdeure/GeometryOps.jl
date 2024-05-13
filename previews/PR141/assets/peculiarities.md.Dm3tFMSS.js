import{_ as e,c as o,o as t,a6 as a}from"./chunks/framework.Ckah3iQP.js";const m=JSON.parse('{"title":"Peculiarities","description":"","frontmatter":{},"headers":[],"relativePath":"peculiarities.md","filePath":"peculiarities.md","lastUpdated":null}'),r={name:"peculiarities.md"},i=a('<h1 id="Peculiarities" tabindex="-1">Peculiarities <a class="header-anchor" href="#Peculiarities" aria-label="Permalink to &quot;Peculiarities {#Peculiarities}&quot;">​</a></h1><h2 id="What-does-apply-return-and-why?" tabindex="-1">What does <code>apply</code> return and why? <a class="header-anchor" href="#What-does-apply-return-and-why?" aria-label="Permalink to &quot;What does `apply` return and why? {#What-does-apply-return-and-why?}&quot;">​</a></h2><p><code>apply</code> returns the target geometries returned by <code>f</code>, whatever type/package they are from, but geometries, features or feature collections that wrapped the target are replaced with GeoInterace.jl wrappers with matching <code>GeoInterface.trait</code> to the originals. All non-geointerface iterables become <code>Array</code>s. Tables.jl compatible tables are converted either back to the original type if a <code>Tables.materializer</code> is defined, and if not then returned as generic <code>NamedTuple</code> column tables (i.e., a NamedTuple of vectors).</p><p>It is recommended for consistency that <code>f</code> returns GeoInterface geometries unless there is a performance/conversion overhead to doing that.</p><h2 id="Why-do-you-want-me-to-provide-a-target-in-set-operations?" tabindex="-1">Why do you want me to provide a <code>target</code> in set operations? <a class="header-anchor" href="#Why-do-you-want-me-to-provide-a-target-in-set-operations?" aria-label="Permalink to &quot;Why do you want me to provide a `target` in set operations? {#Why-do-you-want-me-to-provide-a-target-in-set-operations?}&quot;">​</a></h2><p>In polygon set operations like <code>intersection</code>, <code>difference</code>, and <code>union</code>, many different geometry types may be obtained - depending on the relationship between the polygons. For example, when performing an union on two nonintersecting polygons, one would technically have two disjoint polygons as an output.</p><p>We use the <code>target</code> keyword to allow the user to control which kinds of geometry they want back. For example, setting <code>target</code> to <code>PolygonTrait</code> will cause a vector of polygons to be returned (this is the only currently supported behaviour). In future, we may implement <code>MultiPolygonTrait</code> or <code>GeometryCollectionTrait</code> targets which will return a single geometry, as LibGEOS and ArchGDAL do.</p><p>This also allows for a lot more type stability - when you ask for polygons, we won&#39;t return a geometrycollection with line segments. Especially in simulation workflows, this is excellent for simplified data processing.</p><h2 id="_True-and-_False-(or-BoolsAsTypes)" tabindex="-1"><code>_True</code> and <code>_False</code> (or <code>BoolsAsTypes</code>) <a class="header-anchor" href="#_True-and-_False-(or-BoolsAsTypes)" aria-label="Permalink to &quot;`_True` and `_False` (or `BoolsAsTypes`) {#_True-and-_False-(or-BoolsAsTypes)}&quot;">​</a></h2><div class="warning custom-block"><p class="custom-block-title">Warning</p><p>These are internals and explicitly <em>not</em> public API, meaning they may change at any time!</p></div><p>When dispatch can be controlled by the value of a boolean variable, this introduces type instability. Instead of introducing type instability, we chose to encode our boolean decision variables, like <code>threaded</code> and <code>calc_extent</code> in <code>apply</code>, as types. This allows the compiler to reason about what will happen, and call the correct compiled method, in a stable way without worrying about</p>',11),n=[i];function s(l,c,d,p,h,u){return t(),o("div",null,n)}const g=e(r,[["render",s]]);export{m as __pageData,g as default};
