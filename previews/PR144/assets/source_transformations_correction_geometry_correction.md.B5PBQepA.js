import{_ as i,c as s,o as e,a6 as t}from"./chunks/framework.BZUDmi18.js";const g=JSON.parse('{"title":"Geometry Corrections","description":"","frontmatter":{},"headers":[],"relativePath":"source/transformations/correction/geometry_correction.md","filePath":"source/transformations/correction/geometry_correction.md","lastUpdated":null}'),a={name:"source/transformations/correction/geometry_correction.md"},n=t(`<h1 id="Geometry-Corrections" tabindex="-1">Geometry Corrections <a class="header-anchor" href="#Geometry-Corrections" aria-label="Permalink to &quot;Geometry Corrections {#Geometry-Corrections}&quot;">​</a></h1><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">export</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> fix</span></span></code></pre></div><p>This file simply defines the <code>GeometryCorrection</code> abstract type, and the interface that any <code>GeometryCorrection</code> must implement.</p><p>A geometry correction is a transformation that is applied to a geometry to correct it in some way.</p><p>For example, a <code>ClosedRing</code> correction might be applied to a <code>Polygon</code> to ensure that its exterior ring is closed.</p><h2 id="Interface" tabindex="-1">Interface <a class="header-anchor" href="#Interface" aria-label="Permalink to &quot;Interface {#Interface}&quot;">​</a></h2><p>All <code>GeometryCorrection</code>s are callable structs which, when called, apply the correction to the given geometry, and return either a copy or the original geometry (if nothing needed to be corrected).</p><p>See below for the full interface specification.</p><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="GeometryOps.GeometryCorrection-source-transformations-correction-geometry_correction" href="#GeometryOps.GeometryCorrection-source-transformations-correction-geometry_correction">#</a> <b><u>GeometryOps.GeometryCorrection</u></b> — <i>Type</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">abstract type</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GeometryCorrection</span></span></code></pre></div><p>This abstract type represents a geometry correction.</p><p><strong>Interface</strong></p><p>Any <code>GeometryCorrection</code> must implement two functions: * <code>application_level(::GeometryCorrection)::AbstractGeometryTrait</code>: This function should return the <code>GeoInterface</code> trait that the correction is intended to be applied to, like <code>PointTrait</code> or <code>LineStringTrait</code> or <code>PolygonTrait</code>. * <code>(::GeometryCorrection)(::AbstractGeometryTrait, geometry)::(some_geometry)</code>: This function should apply the correction to the given geometry, and return a new geometry.</p><p><a href="https://github.com/JuliaGeo/GeometryOps.jl/blob/293a403e8ce2120666b26e43fc3ae4dcbac07cfa/src/transformations/correction/geometry_correction.jl#L28-L38" target="_blank" rel="noreferrer">source</a></p></div><br><p>Any geometry correction must implement the interface as given above.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    abstract type GeometryCorrection</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">This abstract type represents a geometry correction.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"># Interface</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Any \`GeometryCorrection\` must implement two functions:</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    * \`application_level(::GeometryCorrection)::AbstractGeometryTrait\`: This function should return the \`GeoInterface\` trait that the correction is intended to be applied to, like \`PointTrait\` or \`LineStringTrait\` or \`PolygonTrait\`.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    * \`(::GeometryCorrection)(::AbstractGeometryTrait, geometry)::(some_geometry)\`: This function should apply the correction to the given geometry, and return a new geometry.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">abstract type</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GeometryCorrection </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">application_level</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(gc</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">GeometryCorrection</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> error</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Not implemented yet for </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$(gc)</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(gc</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">GeometryCorrection</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)(geometry) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> gc</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">trait</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(geometry), geometry)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(gc</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">GeometryCorrection</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)(trait</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">GI.AbstractGeometryTrait</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, geometry) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> error</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Not implemented yet for </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$(gc)</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> and </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$(trait)</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">.&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> fix</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(geometry; corrections </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GeometryCorrection[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ClosedRing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(),], kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    traits </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> application_level</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(corrections)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    final_geometry </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> geometry</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Trait </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">PointTrait, GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">MultiPointTrait, GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">LineStringTrait, GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">LinearRingTrait, GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">MultiLineStringTrait, GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">PolygonTrait, GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">MultiPolygonTrait)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        available_corrections </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> findall</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Trait, traits)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        isempty</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(available_corrections) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;&amp;</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> continue</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        @debug</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;Correcting for </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$(Trait)</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        net_function </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> reduce</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">∘</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, corrections[available_corrections])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        final_geometry </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> apply</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(net_function, Trait, final_geometry; kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    return</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> final_geometry</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><h2 id="Available-corrections" tabindex="-1">Available corrections <a class="header-anchor" href="#Available-corrections" aria-label="Permalink to &quot;Available corrections {#Available-corrections}&quot;">​</a></h2><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="GeometryOps.ClosedRing-source-transformations-correction-geometry_correction" href="#GeometryOps.ClosedRing-source-transformations-correction-geometry_correction">#</a> <b><u>GeometryOps.ClosedRing</u></b> — <i>Type</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ClosedRing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">() </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> GeometryCorrection</span></span></code></pre></div><p>This correction ensures that a polygon&#39;s exterior and interior rings are closed.</p><p>It can be called on any geometry correction as usual.</p><p>See also <a href="/GeometryOps.jl/previews/PR144/api#GeometryOps.GeometryCorrection"><code>GeometryCorrection</code></a>.</p><p><a href="https://github.com/JuliaGeo/GeometryOps.jl/blob/293a403e8ce2120666b26e43fc3ae4dcbac07cfa/src/transformations/correction/closed_ring.jl#L38-L46" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="GeometryOps.DiffIntersectingPolygons-source-transformations-correction-geometry_correction" href="#GeometryOps.DiffIntersectingPolygons-source-transformations-correction-geometry_correction">#</a> <b><u>GeometryOps.DiffIntersectingPolygons</u></b> — <i>Type</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">DiffIntersectingPolygons</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">() </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> GeometryCorrection</span></span></code></pre></div><p>This correction ensures that the polygons included in a multipolygon aren&#39;t intersecting. If any polygon&#39;s are intersecting, they will be made nonintersecting through the <a href="/GeometryOps.jl/previews/PR144/api#GeometryOps.difference-Union{Tuple{T}, Tuple{Any, Any}, Tuple{Any, Any, Type{T}}} where T&lt;:AbstractFloat"><code>difference</code></a> operation to create a unique set of disjoint (other than potentially connections by a single point) polygons covering the same area. See also <a href="/GeometryOps.jl/previews/PR144/api#GeometryOps.GeometryCorrection"><code>GeometryCorrection</code></a>, <a href="/GeometryOps.jl/previews/PR144/api#GeometryOps.UnionIntersectingPolygons"><code>UnionIntersectingPolygons</code></a>.</p><p><a href="https://github.com/JuliaGeo/GeometryOps.jl/blob/293a403e8ce2120666b26e43fc3ae4dcbac07cfa/src/transformations/correction/intersecting_polygons.jl#L92-L99" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="GeometryOps.GeometryCorrection-source-transformations-correction-geometry_correction-2" href="#GeometryOps.GeometryCorrection-source-transformations-correction-geometry_correction-2">#</a> <b><u>GeometryOps.GeometryCorrection</u></b> — <i>Type</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">abstract type</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GeometryCorrection</span></span></code></pre></div><p>This abstract type represents a geometry correction.</p><p><strong>Interface</strong></p><p>Any <code>GeometryCorrection</code> must implement two functions: * <code>application_level(::GeometryCorrection)::AbstractGeometryTrait</code>: This function should return the <code>GeoInterface</code> trait that the correction is intended to be applied to, like <code>PointTrait</code> or <code>LineStringTrait</code> or <code>PolygonTrait</code>. * <code>(::GeometryCorrection)(::AbstractGeometryTrait, geometry)::(some_geometry)</code>: This function should apply the correction to the given geometry, and return a new geometry.</p><p><a href="https://github.com/JuliaGeo/GeometryOps.jl/blob/293a403e8ce2120666b26e43fc3ae4dcbac07cfa/src/transformations/correction/geometry_correction.jl#L28-L38" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="GeometryOps.UnionIntersectingPolygons-source-transformations-correction-geometry_correction" href="#GeometryOps.UnionIntersectingPolygons-source-transformations-correction-geometry_correction">#</a> <b><u>GeometryOps.UnionIntersectingPolygons</u></b> — <i>Type</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">UnionIntersectingPolygons</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">() </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> GeometryCorrection</span></span></code></pre></div><p>This correction ensures that the polygon&#39;s included in a multipolygon aren&#39;t intersecting. If any polygon&#39;s are intersecting, they will be combined through the union operation to create a unique set of disjoint (other than potentially connections by a single point) polygons covering the same area.</p><p>See also <a href="/GeometryOps.jl/previews/PR144/api#GeometryOps.GeometryCorrection"><code>GeometryCorrection</code></a>.</p><p><a href="https://github.com/JuliaGeo/GeometryOps.jl/blob/293a403e8ce2120666b26e43fc3ae4dcbac07cfa/src/transformations/correction/intersecting_polygons.jl#L47-L56" target="_blank" rel="noreferrer">source</a></p></div><br><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,23),r=[n];function o(l,p,h,c,k,d){return e(),s("div",null,r)}const m=i(a,[["render",o]]);export{g as __pageData,m as default};
