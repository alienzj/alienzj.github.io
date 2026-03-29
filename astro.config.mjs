// @ts-check

import mdx from '@astrojs/mdx';
import sitemap from '@astrojs/sitemap';
import { defineConfig } from 'astro/config';

// https://astro.build/config
export default defineConfig({
        site: 'https://alienzj.org',
        integrations: [mdx(), sitemap()],
        markdown: {
                shikiConfig: {
                        themes: {
                                light: 'github-light',
                                dark: 'github-dark',
                        },
                },
        },
        image: {
                domains: ['alienzj.org'],
        },
});

