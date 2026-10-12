// properties.ts — state for the properties sidebar.
//
// Two stores, and they are here rather than in theme.ts or status.ts because
// neither is about a theme or a measurement: one is whether the panel is open,
// the other is a display preference the panel owns.
//
// Plain writables with no persistence, matching darkMode in theme.ts. Nothing in
// this app persists UI state yet, and adding localStorage for one preference
// would make this the only setting that survives a restart -- a surprise, not a
// feature. When persistence arrives it should arrive for all of them at once.
//
// File > Open Recent (src-tauri/src/recent.rs) does persist, and is not the
// exception that breaks this: a recent-files list is document history, not a
// preference -- an Open Recent that forgets on restart is not the feature at all.
// No UI setting has followed it, and none should until they all can.

import { writable } from 'svelte/store';

/** Whether the properties sidebar is showing. Toolbar's Sidebar group toggles it. */
export const propertiesOpen = writable(false);

/** Whether a code cell shows its `In[n]` label once it has been evaluated.
 *
 *  Lives here rather than in the toolbar because it is a persistent preference
 *  about how the notebook reads, not an action -- the toolbar is for verbs. */
export const showExecLabels = writable(true);

/** Global UI scale, applied by App.svelte as the root font size.
 *
 *  Lives here because the properties panel offers it, but it is NOT new: the
 *  Cmd+= / Cmd+- / Cmd+0 bindings have driven this value since the toolbar
 *  landed, as a local in App.svelte. Lifting it to a store is what lets the panel
 *  show the same number the keyboard changes -- a second, independent scale would
 *  drift from the first the moment either was used.
 *
 *  The steps the panel offers are exactly representable in binary (0.75, 1, 1.25,
 *  1.5), so a comparison against one is an exact equality rather than an epsilon.
 *  The keyboard still moves in 0.1 increments, which is why the panel highlights a
 *  step only when the value is exactly on it. */
export const uiScale = writable(1.0);

/** The steps the panel offers. Not the keyboard's 0.1: a discrete set is what a
 *  set of buttons can honestly represent. */
export const UI_SCALE_STEPS = [0.75, 1, 1.25, 1.5];

/** Whether code cells offer the builtin/session-symbol completion dropdown.
 *  Read by the completion source and the Tab binding in lib/completion.ts; the
 *  Properties panel toggles it. Plain writable, no persistence — same decision
 *  as the preferences above. */
export const autocompleteEnabled = writable(true);

/** The active syntax-highlighting colour scheme id (see lib/schemes.ts).
 *  'default' is the adaptive Catppuccin light/dark palette; the others are fixed
 *  palettes applied over it. App.svelte applies this; the panel chooses it. */
export const colorScheme = writable('default');

/** View › Zoom: the magnification of every cell's TEXT -- inputs, outputs, text and heading
 *  cells, In/Out labels -- across all notebooks, as a factor (1 = 100%). Distinct from
 *  `uiScale`, which scales the whole window (toolbar and panels too): this one leaves the
 *  chrome alone, like a notebook's magnification. App.svelte publishes it as the CSS variable
 *  --cell-zoom, which the cell styles multiply their font sizes by. Plain writable, no
 *  persistence — same decision as the preferences above. */
export const cellZoom = writable(1.0);

/** The Zoom menu's choices: 50% to 200% in 10% steps. Integers (percent) so a menu id and a
 *  comparison are exact; the factor is percent / 100. */
export const CELL_ZOOM_PERCENTS = [50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200];

/** The notation every expression result is displayed in (View › Output Format):
 *  'standard' — StandardForm, Mathilda's own names typeset, the default;
 *  'traditional' — TraditionalForm, textbook notation;
 *  'input' — InputForm, 1D re-readable text.
 *  A result switched with its own right-click Convert To keeps that choice; this is the
 *  form of every result that has not been. Plain writable, no persistence — same decision
 *  as the preferences above. */
export type OutputFormat = 'input' | 'standard' | 'traditional';
export const outputFormat = writable<OutputFormat>('standard');

