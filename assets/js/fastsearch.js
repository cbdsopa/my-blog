import * as params from '@params';

const resList = document.getElementById('searchResults');
const sInput = document.getElementById('searchInput');
const searchBox = document.getElementById('searchbox');

let fuse;
let currentElement = null;
let firstResult = null;
let lastResult = null;

const defaultFuseOptions = {
    distance: 100,
    threshold: 0.4,
    ignoreLocation: true,
    includeMatches: true,
    minMatchCharLength: 2,
    keys: ['title', 'permalink', 'summary', 'content']
};

const buildFuseOptions = () => {
    if (!params.fuseOpts) {
        return defaultFuseOptions;
    }

    return {
        isCaseSensitive: params.fuseOpts.iscasesensitive ?? false,
        includeScore: params.fuseOpts.includescore ?? false,
        includeMatches: true,
        minMatchCharLength: params.fuseOpts.minmatchcharlength ?? 2,
        shouldSort: params.fuseOpts.shouldsort ?? true,
        findAllMatches: params.fuseOpts.findallmatches ?? true,
        keys: params.fuseOpts.keys ?? defaultFuseOptions.keys,
        location: params.fuseOpts.location ?? 0,
        threshold: params.fuseOpts.threshold ?? defaultFuseOptions.threshold,
        distance: params.fuseOpts.distance ?? defaultFuseOptions.distance,
        ignoreLocation: params.fuseOpts.ignorelocation ?? defaultFuseOptions.ignoreLocation
    };
};

const debounce = (fn, delay) => {
    let timeout;
    return (...args) => {
        clearTimeout(timeout);
        timeout = window.setTimeout(() => fn(...args), delay);
    };
};

const reset = () => {
    currentElement = null;
    firstResult = null;
    lastResult = null;
    resList.innerHTML = '';
    sInput.value = '';
    sInput.focus();
};

const setActiveResult = (element) => {
    document.querySelectorAll('.focus').forEach((item) => item.classList.remove('focus'));

    if (!element) {
        return;
    }

    element.focus();
    element.parentElement?.classList.add('focus');
    currentElement = element;
};

// Extract a snippet from content around matched positions.
// Returns { text, adjustedIndices } where adjustedIndices are relative to the snippet.
const getSnippet = (content, matches, maxLen) => {
    maxLen = maxLen || 120;
    const empty = { text: '', adjustedIndices: [] };

    if (!content || !matches || matches.length === 0) {
        return empty;
    }

    // Gather all match indices across all keys, filtering out non-content matches
    const contentMatch = matches.find(m => m.key === 'content' && m.indices && m.indices.length > 0);
    const summaryMatch = matches.find(m => m.key === 'summary' && m.indices && m.indices.length > 0);
    const titleMatch = matches.find(m => m.key === 'title' && m.indices && m.indices.length > 0);

    const targetMatch = contentMatch || summaryMatch || titleMatch;
    if (!targetMatch) {
        return { text: content.substring(0, maxLen) + (content.length > maxLen ? '…' : ''), adjustedIndices: [] };
    }

    // Use the first index from the best match as the anchor
    const idx = targetMatch.indices[0];
    const matchStart = idx[0];
    const matchEnd = idx[1];
    const matchLen = matchEnd - matchStart;

    // Expand context around the match
    const context = Math.max(20, Math.floor((maxLen - matchLen) / 2));
    const start = Math.max(0, matchStart - context);
    const end = Math.min(content.length, matchEnd + context);

    let snippet = content.substring(start, end);
    if (start > 0) snippet = '…' + snippet;
    if (end < content.length) snippet = snippet + '…';

    // Adjust match indices relative to the snippet
    const ellipsisOffset = start > 0 ? 1 : 0; // '…' is one char
    const adjustedIndices = [];

    // Collect all match indices, adjusting them to snippet coordinates
    for (const match of matches) {
        if (!match.indices) continue;
        for (const [ms, me] of match.indices) {
            // Only include if this match overlaps with our snippet range
            if (me <= start || ms >= end) continue;
            const adjStart = Math.max(0, ms - start) + ellipsisOffset;
            const adjEnd = Math.min(snippet.length, me - start) + ellipsisOffset;
            if (adjStart < adjEnd && adjStart < snippet.length) {
                adjustedIndices.push([adjStart, adjEnd]);
            }
        }
    }

    return { text: snippet, adjustedIndices };
};

// Highlight matched terms in a snippet using adjusted indices
const highlightMatches = (text, adjustedIndices) => {
    if (!text || !adjustedIndices || adjustedIndices.length === 0) {
        return text;
    }

    // Build a character-level set of positions to highlight
    const positions = new Set();
    for (const [s, e] of adjustedIndices) {
        for (let i = s; i < e && i < text.length; i++) {
            positions.add(i);
        }
    }

    if (positions.size === 0) return text;

    // Build highlighted HTML
    let result = '';
    let inMark = false;
    for (let i = 0; i < text.length; i++) {
        if (positions.has(i)) {
            if (!inMark) {
                result += '<mark>';
                inMark = true;
            }
        } else {
            if (inMark) {
                result += '</mark>';
                inMark = false;
            }
        }
        result += text[i];
    }
    if (inMark) result += '</mark>';

    return result;
};

const renderResults = (results) => {
    if (!Array.isArray(results) || results.length === 0) {
        resList.innerHTML = '';
        firstResult = lastResult = currentElement = null;
        return;
    }

    const fragment = document.createDocumentFragment();

    for (const result of results) {
        const li = document.createElement('li');
        li.className = 'search-result-item';

        // Title line
        const titleDiv = document.createElement('div');
        titleDiv.className = 'search-result-title';
        titleDiv.textContent = result.item.title;

        // Preview snippet with highlighted matches
        const previewDiv = document.createElement('div');
        previewDiv.className = 'search-result-preview';

        const content = result.item.content || '';
        const { text: snippet, adjustedIndices } = getSnippet(content, result.matches, 150);
        const highlighted = highlightMatches(snippet, adjustedIndices);
        previewDiv.innerHTML = highlighted;

        // Chevron icon
        const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
        svg.setAttribute('width', '24');
        svg.setAttribute('height', '24');
        svg.setAttribute('viewBox', '0 0 24 24');
        svg.setAttribute('fill', 'none');
        svg.setAttribute('stroke', 'currentColor');
        svg.setAttribute('stroke-width', '2');
        svg.setAttribute('stroke-linecap', 'round');
        svg.setAttribute('stroke-linejoin', 'round');
        svg.classList.add('feather', 'feather-chevrons-right');
        svg.innerHTML = '<polyline points="13 17 18 12 13 7"></polyline><polyline points="6 17 11 12 6 7"></polyline>';

        const link = document.createElement('a');
        link.className = 'entry-link';
        link.href = result.item.permalink;
        link.setAttribute('aria-label', result.item.title);

        li.appendChild(titleDiv);
        li.appendChild(previewDiv);
        li.appendChild(svg);
        li.appendChild(link);
        fragment.appendChild(li);
    }

    resList.innerHTML = '';
    resList.appendChild(fragment);
    firstResult = resList.firstElementChild;
    lastResult = resList.lastElementChild;
};

const performSearch = () => {
    if (!fuse) {
        return;
    }

    const query = sInput.value.trim();
    if (!query) {
        renderResults([]);
        return;
    }

    const searchOptions = params.fuseOpts?.limit ? { limit: params.fuseOpts.limit } : undefined;
    const results = searchOptions ? fuse.search(query, searchOptions) : fuse.search(query);
    renderResults(results);
};

const initSearch = async () => {
    if (!sInput || !resList) {
        return;
    }

    sInput.disabled = false;
    sInput.focus();

    try {
        const response = await fetch('../index.json');
        if (!response.ok) {
            throw new Error(`Search index load failed: ${response.status}`);
        }

        const data = await response.json();
        if (data) {
            fuse = new Fuse(data, buildFuseOptions());
        }
    } catch (error) {
        console.error(error);
    }
};

window.addEventListener('load', initSearch);

sInput?.addEventListener('input', debounce(performSearch, 150));

sInput?.addEventListener('search', () => {
    if (!sInput.value) {
        reset();
    }
});

document.addEventListener('keydown', (event) => {
    const { key } = event;
    const active = document.activeElement;
    const isInSearchBox = searchBox?.contains(active);

    if (key === 'Escape') {
        reset();
        return;
    }

    if (!firstResult || !isInSearchBox) {
        return;
    }

    if (key === 'ArrowDown') {
        event.preventDefault();

        if (active === sInput) {
            setActiveResult(firstResult.querySelector('.entry-link'));
        } else if (active?.parentElement !== lastResult) {
            setActiveResult(active?.parentElement?.nextElementSibling?.querySelector('.entry-link'));
        }
    } else if (key === 'ArrowUp') {
        event.preventDefault();

        if (active?.parentElement === firstResult) {
            setActiveResult(sInput);
        } else if (active !== sInput) {
            setActiveResult(active?.parentElement?.previousElementSibling?.querySelector('.entry-link'));
        }
    } else if (key === 'ArrowRight') {
        if (active?.matches?.('.entry-link')) {
            active.click();
        }
    }
});
