import{U as o,a8 as p,a9 as u,aa as l,ab as c,ac as f,ad as d,ae as m,af as h,ag as g,ah as A,d as P,u as v,y,x as w,ai as C,aj as R,ak as b,a6 as E}from"./chunks/framework.yROC2_57.js";import{R as S}from"./chunks/theme.DxMxEAIk.js";function i(e){if(e.extends){const a=i(e.extends);return{...a,...e,async enhanceApp(t){a.enhanceApp&&await a.enhanceApp(t),e.enhanceApp&&await e.enhanceApp(t)}}}return e}const s=i(S),T=P({name:"VitePressApp",setup(){const{site:e,lang:a,dir:t}=v();return y(()=>{w(()=>{document.documentElement.lang=a.value,document.documentElement.dir=t.value})}),e.value.router.prefetchLinks&&C(),R(),b(),s.setup&&s.setup(),()=>E(s.Layout)}});async function _(){globalThis.__VITEPRESS__=!0;const e=x(),a=j();a.provide(u,e);const t=l(e.route);return a.provide(c,t),a.component("Content",f),a.component("ClientOnly",d),Object.defineProperties(a.config.globalProperties,{$frontmatter:{get(){return t.frontmatter.value}},$params:{get(){return t.page.value.params}}}),s.enhanceApp&&await s.enhanceApp({app:a,router:e,siteData:m}),{app:a,router:e,data:t}}function j(){return h(T)}function x(){let e=o,a;return g(t=>{let n=A(t),r=null;return n&&(e&&(a=n),(e||a===n)&&(n=n.replace(/\.js$/,".lean.js")),r=import(n)),o&&(e=!1),r},s.NotFound)}o&&_().then(({app:e,router:a,data:t})=>{a.go().then(()=>{p(a.route,t.site),e.mount("#app")})});export{_ as createApp};
