## Gufo

tools to fit multiple gaussians in datacubes

### Installation

Add `<path to radiobs>` to `PYTHONPATH`
  
  - in your `.cshrc`
  `setenv PYTHONPATH <path_to_gufo>:${PYTHONPATH}`
  - in your `.bashrc`
  `export PYTHONPATH=<path_to_gufo>:$PYTHONPATH`
  
### Usage
 In `python` environment type `import gufo`
 
 find which classes are available.
 
 `gufo.__dict__.keys()`
 
 We choose to load the table containind the information about the lines to:
 
 `tp= gufo.tPlay()`
 
  use the module:
 `lineInfo = tp.openLineList()`

 
