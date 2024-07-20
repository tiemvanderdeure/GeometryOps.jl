import{_ as s,c as i,o as a,a7 as e}from"./chunks/framework.BcK3SVj7.js";const y=JSON.parse('{"title":"Geometry reprojection","description":"","frontmatter":{},"headers":[],"relativePath":"source/transformations/reproject.md","filePath":"source/transformations/reproject.md","lastUpdated":null}'),n={name:"source/transformations/reproject.md"},t=e('<h1 id="Geometry-reprojection" tabindex="-1">Geometry reprojection <a class="header-anchor" href="#Geometry-reprojection" aria-label="Permalink to &quot;Geometry reprojection {#Geometry-reprojection}&quot;">​</a></h1><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">export</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> reproject</span></span></code></pre></div><p>This file is pretty simple - it simply reprojects a geometry pointwise from one CRS to another. It uses the <code>Proj</code> package for the transformation, but this could be moved to an extension if needed.</p><p>Note that the actual implementation is in the <code>GeometryOpsProjExt</code> extension module.</p><p>This works using the <code>apply</code> functionality.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    reproject(geometry; source_crs, target_crs, transform, always_xy, time)</span></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    reproject(geometry, source_crs, target_crs; always_xy, time)</span></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    reproject(geometry, transform; always_xy, time)</span></span>\n<span class="line"></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Reproject any GeoInterface.jl compatible `geometry` from `source_crs` to `target_crs`.</span></span>\n<span class="line"></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">The returned object will be constructed from `GeoInterface.WrapperGeometry`</span></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">geometries, wrapping views of a `Vector{Proj.Point{D}}`, where `D` is the dimension.</span></span>\n<span class="line"></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">!!! tip</span></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    The `Proj.jl` package must be loaded for this method to work,</span></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    since it is implemented in a package extension.</span></span>\n<span class="line"></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"># Arguments</span></span>\n<span class="line"></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">- `geometry`: Any GeoInterface.jl compatible geometries.</span></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">- `source_crs`: the source coordinate reference system, as a GeoFormatTypes.jl object or a string.</span></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">- `target_crs`: the target coordinate reference system, as a GeoFormatTypes.jl object or a string.</span></span>\n<span class="line"></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">If these a passed as keywords, `transform` will take priority.</span></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Without it `target_crs` is always needed, and `source_crs` is</span></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">needed if it is not retrievable from the geometry with `GeoInterface.crs(geometry)`.</span></span>\n<span class="line"></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"># Keywords</span></span>\n<span class="line"></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">- `always_xy`: force x, y coordinate order, `true` by default.</span></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    `false` will expect and return points in the crs coordinate order.</span></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">- `time`: the time for the coordinates. `Inf` by default.</span></span>\n<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$APPLY_KEYWORDS</span></span>\n<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>\n<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> reproject </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><h2 id="Method-error-handling" tabindex="-1">Method error handling <a class="header-anchor" href="#Method-error-handling" aria-label="Permalink to &quot;Method error handling {#Method-error-handling}&quot;">​</a></h2><p>We also inject a method error handler, which prints a suggestion if the Proj extension is not loaded.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> _reproject_error_hinter</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, exc, argtypes, kwargs)</span></span>\n<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> isnothing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(Base</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">get_extension</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(GeometryOps, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:GeometryOpsProjExt</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> exc</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">f </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> reproject</span></span>\n<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        print</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\n\\n</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">The `reproject` method requires the Proj.jl package to be explicitly loaded.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\n</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>\n<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        print</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;You can do this by simply typing &quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>\n<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        printstyled</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;using Proj&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :cyan</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, bold </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>\n<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        println</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot; in your REPL, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\n</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">or otherwise loading Proj.jl via using or import.&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>\n<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    else</span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"> # this is a more general error</span></span>\n<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        nothing</span></span>\n<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>\n<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>',11),p=[t];function l(r,h,o,k,d,c){return a(),i("div",null,p)}const F=s(n,[["render",l]]);export{y as __pageData,F as default};
