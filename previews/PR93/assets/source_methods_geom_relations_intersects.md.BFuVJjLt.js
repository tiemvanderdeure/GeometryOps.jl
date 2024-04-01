import{_ as s,c as i,o as a,a7 as n}from"./chunks/framework.DYf-6cjo.js";const e="/GeometryOps.jl/previews/PR93/assets/zznxfpz.BwVBXAhv.png",y=JSON.parse('{"title":"Intersection checks","description":"","frontmatter":{},"headers":[],"relativePath":"source/methods/geom_relations/intersects.md","filePath":"source/methods/geom_relations/intersects.md","lastUpdated":null}'),t={name:"source/methods/geom_relations/intersects.md"},l=n(`<h1 id="Intersection-checks" tabindex="-1">Intersection checks <a class="header-anchor" href="#Intersection-checks" aria-label="Permalink to &quot;Intersection checks {#Intersection-checks}&quot;">​</a></h1><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">export</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> intersects</span></span></code></pre></div><h2 id="What-is-intersects?" tabindex="-1">What is <code>intersects</code>? <a class="header-anchor" href="#What-is-intersects?" aria-label="Permalink to &quot;What is \`intersects\`? {#What-is-intersects?}&quot;">​</a></h2><p>The intersects function checks if a given geometry intersects with another geometry, or in other words, the either the interiors or boundaries of the two geometries intersect.</p><p>To provide an example, consider these two lines:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">import</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GeometryOps </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">as</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GO</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">import</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GeoInterface </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">as</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GI</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Makie</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> CairoMakie</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">line1 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Line</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">([(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">124.584961</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">12.768946</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">126.738281</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">17.224758</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">line2 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Line</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">([(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">123.354492</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">15.961329</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">127.22168</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">14.008696</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">f, a, p </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> lines</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">getpoint</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(line1))</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">lines!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">getpoint</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(line2))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">f</span></span></code></pre></div><p><img src="`+e+`" alt=""></p><p>We can see that they intersect, so we expect intersects to return true, and we can visualize the intersection point in red.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">GO</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">intersects</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(line1, line2)  </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># true</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span>true</span></span></code></pre></div><h2 id="Implementation" tabindex="-1">Implementation <a class="header-anchor" href="#Implementation" aria-label="Permalink to &quot;Implementation {#Implementation}&quot;">​</a></h2><p>This is the GeoInterface-compatible implementation.</p><p>Given that intersects is the exact opposite of disjoint, we simply pass the two inputs variables, swapped in order, to disjoint.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    intersects(geom1, geom2)::Bool</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Return true if the interiors or boundaries of the two geometries interact.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\`intersects\` returns the exact opposite result of \`disjoint\`.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"># Example</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\`\`\`jldoctest</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">import GeoInterface as GI, GeometryOps as GO</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">line1 = GI.Line([(124.584961,-12.768946), (126.738281,-17.224758)])</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">line2 = GI.Line([(123.354492,-15.961329), (127.22168,-14.008696)])</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">GO.intersects(line1, line2)</span></span></code></pre></div><p>output</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">true</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\`\`\`</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">intersects(geom1, geom2) = !disjoint(geom1, geom2)</span></span>
<span class="line"></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">#= Returns true if there is at least one intersection between edges within the</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">two lists of edges. =#</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">function _line_intersects(</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    edges_a::Vector{&lt;:Edge},</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    edges_b::Vector{&lt;:Edge}</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">)</span></span></code></pre></div><p>Extents.intersects(to_extent(edges_a), to_extent(edges_b)) || return false</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> edge_a </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> edges_a</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> edge_b </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> edges_b</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">            _line_intersects</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(edge_a, edge_b) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;&amp;</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> return</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    return</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> false</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><p>Returns true if there is at least one intersection between two edges.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> _line_intersects</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(edge_a</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Edge</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, edge_b</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Edge</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    seg_val, _, _ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> _intersection_point</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(Float64, edge_a, edge_b)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    return</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> seg_val </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">!=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> line_out</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,22),p=[l];function h(k,r,d,g,o,c){return a(),i("div",null,p)}const F=s(t,[["render",h]]);export{y as __pageData,F as default};
