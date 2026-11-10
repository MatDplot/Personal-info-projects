# Projets-perso commands (bash and zsh), sourced by Setup's dotfiles/shell_common.sh; listed in COMMANDS.md, shown by Mhelp_site.
PROJETS_REPO=${PROJETS_REPO:-$HOME/Code/Perso/Projets-perso}

# Each project's scripts from its own folder, one at a time (the pgf figures are heavy on RAM).
# pmssm-sampling/generate.py needs SOFTSUSY and molecular-dynamics needs its C++ results: left out.
projets_run() {
    local f failed=""
    [ $# -gt 0 ] || set -- $(cd "$PROJETS_REPO" && ls -d */*/ | grep -v '^_' | sed 's:/$::')
    for f in $(cd "$PROJETS_REPO" && find "$@" -maxdepth 1 -name '*.py' ! -name generate.py | sort); do
        printf '\n━━ %s ━━\n' "$f"
        ( cd "$PROJETS_REPO/${f%/*}" && python3 "${f##*/}" ) || failed="$failed $f"
    done
    [ -z "$failed" ] && echo "✓ Tout passe" || { echo "✘ Échec :$failed"; return 1; }
}
