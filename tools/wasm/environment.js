// The environment variables of a test: node's, or those a page or shell
// runner sets in globalThis.FORTSPARSE_ENV.
Module.preRun = (Module.preRun || []).concat([() => {
  if (typeof ENV === "undefined") return;
  if (typeof process !== "undefined" && process.env) Object.assign(ENV, process.env);
  if (typeof globalThis.FORTSPARSE_ENV !== "undefined") Object.assign(ENV, globalThis.FORTSPARSE_ENV);
}]);
