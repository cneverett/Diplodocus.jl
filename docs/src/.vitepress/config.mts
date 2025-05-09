import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import path from 'path'

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: '/cneverett.github.io\DiplodocusDocs/',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: [
{ text: 'Overview', collapsed: false, items: [
{ text: 'Overview', link: '/Overview/overview' },
{ text: 'Installation', link: '/Overview/installation' },
{ text: 'Examples', collapsed: false, items: [
{ text: 'Examples', link: '/Overview/Examples/examples' }]
 }]
 },
{ text: 'Collisions', collapsed: false, items: [
{ text: 'Overview', link: '/DiplodocusCollisions/overview' }]
 },
{ text: 'Transport', collapsed: false, items: [
{ text: 'Overview', link: '/DiplodocusTransport/overview' }]
 },
{ text: 'Plots', collapsed: false, items: [
{ text: 'Overview', link: '/DiplodocusPlots/overview' }]
 }
]
,
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/cneverett.github.io\DiplodocusDocs/',// TODO: replace this in makedocs!
  title: 'DiplodocusDocs',
  description: 'Documentation for DiplodocusDocs',
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../final_site', // This is required for MarkdownVitepress to work correctly...
  head: [
    
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  
  vite: {
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },
  themeConfig: {
    outline: 'deep',
    
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
{ text: 'Overview', collapsed: false, items: [
{ text: 'Overview', link: '/Overview/overview' },
{ text: 'Installation', link: '/Overview/installation' },
{ text: 'Examples', collapsed: false, items: [
{ text: 'Examples', link: '/Overview/Examples/examples' }]
 }]
 },
{ text: 'Collisions', collapsed: false, items: [
{ text: 'Overview', link: '/DiplodocusCollisions/overview' }]
 },
{ text: 'Transport', collapsed: false, items: [
{ text: 'Overview', link: '/DiplodocusTransport/overview' }]
 },
{ text: 'Plots', collapsed: false, items: [
{ text: 'Overview', link: '/DiplodocusPlots/overview' }]
 }
]
,
    editLink: { pattern: "https://https://github.com/cneverett/DiplodocusDocs/edit/main/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/cneverett/DiplodocusDocs' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
