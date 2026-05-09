/** Category config for dotfiles docs sidebar */
export interface DocMeta {
  section: string;
  title: string;
}

const docs: Record<string, DocMeta> = {
  // Root docs
  README: { section: 'Architecture', title: 'README — Project Overview' },
  CLAUDE: { section: 'Architecture', title: 'CLAUDE.md — Agent Instructions' },
  CHANGELOG: { section: 'Architecture', title: 'Changelog' },

  // Architecture
  'nix-expressions': { section: 'Architecture', title: 'Nix Expressions & Module System' },
  packages: { section: 'Architecture', title: 'Package Hierarchy & Derivations' },
  toolchain: { section: 'Architecture', title: 'Hey Toolchain & CLI' },

  // Desktop
  desktop: { section: 'Desktop', title: 'Desktop Environment & Compositors' },
  themes: { section: 'Desktop', title: 'Theme System & Stylix' },
  software: { section: 'Desktop', title: 'Software Jail & Fake Home' },
  keybindings: { section: 'Desktop', title: 'Keybindings & Input' },
  tmux: { section: 'Desktop', title: 'Tmux Configuration' },

  // Editors
  editors: { section: 'Editors', title: 'Multi-Editor Setup' },
  neovim: { section: 'Editors', title: 'Neovim LSP/DAP' },
  'ai-language-idioms': { section: 'Editors', title: 'AI Language Idioms Guide' },
  'ai-parallel-workflow': { section: 'Editors', title: 'AI Parallel Workflow' },

  // Security
  security: { section: 'Security', title: 'Security Overview' },
  'security-hardening': { section: 'Security', title: 'Security Hardening' },
  'security-auth-logic': { section: 'Security', title: 'Authentication Logic' },
  'vulnerability-response': { section: 'Security', title: 'Vulnerability Response' },
  'git-signing': { section: 'Security', title: 'Git SSH Signing' },

  // Networking + SSH
  'networking-vpn': { section: 'Networking', title: 'VPN & Tailscale' },
  'networking-proxy': { section: 'Networking', title: 'Proxy & Sing-box' },
  ssh: { section: 'Networking', title: 'SSH Architecture' },
  'web-services': { section: 'Networking', title: 'Web Services & Nginx' },

  // Services
  'sso-identity': { section: 'Services', title: 'SSO & Identity Management' },
  'containers-virt': { section: 'Services', title: 'Containers & Virtualization' },
  'systemd-services': { section: 'Services', title: 'Systemd Services' },

  // Storage & Hardware
  'storage-disko': { section: 'Storage & Hardware', title: 'Storage & Disko' },
  persistence: { section: 'Storage & Hardware', title: 'Persistence & Impermanence' },
  hardware: { section: 'Storage & Hardware', title: 'Hardware Profiles' },
  'power-management': { section: 'Storage & Hardware', title: 'Power Management' },

  // Data & AI
  'data-science': { section: 'Data & AI', title: 'Data Science Suite' },
  'ai-ml': { section: 'Data & AI', title: 'AI/ML Infrastructure' },

  // Operations
  ops: { section: 'Operations', title: 'Operations Guide' },
  workflow: { section: 'Operations', title: 'Development Workflow' },
  'refactor-plan': { section: 'Operations', title: 'Refactoring Plan' },
  'sbc-opi5p': { section: 'Operations', title: 'SBC Orange Pi 5 Plus' },
};

/** Section display order */
export const sectionOrder = [
  'Architecture',
  'Desktop',
  'Editors',
  'Security',
  'Networking',
  'Services',
  'Storage & Hardware',
  'Data & AI',
  'Hosts',
  'Operations',
];

export function getDocMeta(slug: string): DocMeta {
  if (docs[slug]) return docs[slug];
  // Hosts: slug like "hosts/bio-alpha"
  if (slug.startsWith('hosts/')) {
    const hostname = slug.replace('hosts/', '');
    return { section: 'Hosts', title: `Host: ${hostname}` };
  }
  // Fallback: derive title from filename
  const name = slug.split('/').pop() || slug;
  return { section: 'Other', title: name.replace(/-/g, ' ').replace(/\b\w/g, (c) => c.toUpperCase()) };
}

export function groupBySection(
  entries: { slug: string }[],
): { section: string; items: { slug: string; title: string }[] }[] {
  const map = new Map<string, { slug: string; title: string }[]>();
  for (const e of entries) {
    const meta = getDocMeta(e.slug);
    if (!map.has(meta.section)) map.set(meta.section, []);
    map.get(meta.section)!.push({ slug: e.slug, title: meta.title });
  }
  return sectionOrder
    .filter((sec) => map.has(sec))
    .map((section) => ({ section, items: map.get(section)! }))
    .concat(
      [...map.entries()]
        .filter(([sec]) => !sectionOrder.includes(sec))
        .map(([section, items]) => ({ section, items })),
    );
}
