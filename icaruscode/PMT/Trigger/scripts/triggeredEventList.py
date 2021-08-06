#!/usr/bin/env python3

import time

__doc__ = """Prints the list of events passing the specified trigger."""
__author__ = "Gianluca Petrillo (petrillo@slac.stanford.edu)"
__date__ = time.strptime("July 11, 2021", "%B %d, %Y")
__version__ = "1.0"

import sys
import galleryUtils
import ROOT
import logging

DefaultTriggerModuleLabel = "trigger" # may be overridden by command line

TriggerDataProductClass = 'std::vector<raw::Trigger>'

# ------------------------------------------------------------------------------
class EventTagCache:
  def __init__(self,
               run=False, subRun=False, event=True, sourceFile=False,
               human=False, short=False,
               ):
    self.eventTag = None
    self.run = run
    self.subrun = subRun
    self.event = event
    self.sourceFile = sourceFile
    self.human = human
    self.short = short
  # __init__()
  
  def prepareFor(self, event):
    self.eventObj = event
    self.eventTag = None
  
  def __call__(self, event = None):
    if self.eventTag is None or event is not self.eventObj:
      self.eventTag = self._makeTag(event)
    return self.eventTag
  # __call__()
  
  def _makeTag(self, event=None):
    if event is not None: self.eventObj = event
    assert self.eventObj, "Needs to prepareFor(event) this object first!"
    eventAux = self.eventObj.eventAuxiliary()
    eventID = []
    if self.run:
      run = eventAux.run()
      if self.human: eventID.append(f"run {run}")
      elif self.short: eventID.append(f"R:{run}")
      else: eventID.append(str(run))
    # if print run
    if self.subrun:
      subrun = eventAux.subRun()
      if self.human: eventID.append(f"subrun {subrun}")
      elif self.short: eventID.append(f"S:{subrun}")
      else: eventID.append(str(subrun))
    # if print subrun
    if self.event:
      evt = eventAux.event()
      if self.human: eventID.append(f"event {evt}")
      elif self.short: eventID.append(f"E:{evt}")
      else: eventID.append(str(evt))
    # if print event
    if self.sourceFile:
      sourceFile = self.eventObj.getTFile()
      fileName = sourceFile.GetName() if sourceFile else "<unknown>"
      eventInFile = self.eventObj.eventEntry()
      if self.human: eventID.append(f"from file {fileName} entry {eventInFile}")
      elif self.short: eventID.append(f"F:{fileName!r}@{eventInFile}")
      else:
        eventID.append(repr(fileName))
        eventID.append(str(eventInFile))
    # if print event
    return " ".join(eventID)
  # _makeTag()
  
# class EventTagCache


# ------------------------------------------------------------------------------
def findTriggeredEvents(sampleEvents, triggers, args):
  
  # digest the trigger info options
  trigSpecs = triggers.copy()
  for trigSpec in trigSpecs:
    tag = trigSpec['spec']
    if not isinstance(tag, ROOT.art.InputTag):
      if ':' not in tag: tag = DefaultTriggerModuleLabel + ':' + tag
      tag = ROOT.art.InputTag(tag)
    # if
    trigSpec.update({
      'tag': tag,
      'dest': (open(trigSpec['destList'], 'w')
        if trigSpec['destList'] else sys.stdout),
      'count': 0,
      })
  # for
  
  try: # using a context manager here is more complicate...
    
    getTrigger \
     = galleryUtils.make_getValidHandle(TriggerDataProductClass, sampleEvents)
    
    eventTag = EventTagCache(
      run=args.run,
      subRun=args.subrun,
      event=args.event,
      sourceFile=args.sourcefile,
      human=args.human,
      short=args.short,
      )
    
    for iEvent, event in enumerate(galleryUtils.forEach(sampleEvents)):
      
      if args.maxEvents and iEvent >= args.maxEvents:
        logging.warning(f"Limit of {iEvent} events reached: stopping now.")
        break
      
      eventTag.prepareFor(event)
      
      for trigSpec in trigSpecs:
        triggers = getTrigger(trigSpec['tag']).product()
        if len(triggers) == 0: continue
        print(eventTag(), file=trigSpec['dest'])
        trigSpec['count'] += 1
      # for
    
    # for all events
    else: iEvent += 1
    
  finally:
    for trigSpec in trigSpecs:
      if trigSpec['dest'] is not sys.stdout: trigSpec['dest'].close()
    # for
  # try ... finally
  
  for trigSpec in trigSpecs:
    print(
      f"{trigSpec['tag'].encode()}:"
      f" {trigSpec['count']}/{iEvent} events triggered"
      f" ({trigSpec['count']/iEvent*100:g}%)",
      end="",
      )
    if (trigSpec['dest'] is not sys.stdout):
      print(f" (->'{trigSpec['dest'].name}')", end="")
    print()
  # for
  
  return 0
# findTriggeredEvents()


def getAllTriggerTags(sampleEvents):
  
  return sorted(
    sampleEvents.getInputTags[TriggerDataProductClass](),
    key=ROOT.art.InputTag.encode
    )
  
# getAllTriggerTags()


def listTriggerDataProducts(sampleEvents):
  # sorted() also converts a returned `std::vector` into a Python list
  inputTags = getAllTriggerTags(sampleEvents)
  nProducts = len(inputTags)
  print(f"{nProducts} '{TriggerDataProductClass}' data products found:")
  PadLength=len(str(nProducts))
  for iProduct, inputTag in enumerate(inputTags):
    print(f"[#{iProduct+1:0>{PadLength}}] {inputTag.encode()}")
  # for
  return 0
# listTriggerDataProducts()


def buildDefaultOutputFile(spec):
  if isinstance(spec, ROOT.art.InputTag): spec = spec.encode()
  return f"triggeredEvents_{spec.replace(':', '_')}.log"
# buildDefaultOutputFile()


if __name__ == "__main__":
  
  # --- BEGIN -- argument parsing  ---------------------------------------------
  import argparse
  
  parser = argparse.ArgumentParser(description=__doc__)
  
  # positional
  parser.add_argument("InputFile", help="art/ROOT input file or file list")
  parser.add_argument("TriggerSpec", nargs="*",
    help=f"trigger to check: data product name (`{DefaultTriggerModuleLabel}:inst`)"
     f" or just instance name (`M5S10`, becoming `{DefaultTriggerModuleLabel}:M5S10`)"
     )
  
  # options
  inputGroup = parser.add_argument_group("Input options")
  inputGroup.add_argument("--maxEvents", "-n", type=int,
    help="do not process more that this number of events")
  inputGroup.add_argument("--deflabel", type=str, default=DefaultTriggerModuleLabel,
    help="trigger module label used in TriggerSpec by default [%(default)s]")
  inputGroup.add_argument("--all", action="store_true",
    help="processes all available trigger data products")
  inputGroup.add_argument("--all-to-files", dest="allToFiles",
    action="store_true",
    help="like `--all`, redirecting each trigger into its default-named file"
    )
  
  outputFormatGroup = parser.add_argument_group("Output format")
  outputFormatGroup.add_argument("--run", "-r", action="store_true",
    help="print the run number for triggering events [%(default)s]")
  outputFormatGroup.add_argument("--subrun", "-a", action="store_true",
    help="print the subrun number for triggering events [%(default)s]")
  outputFormatGroup.add_argument("--event", "-e", action="store_false",
    help="print the event number for triggering events [%(default)s]")
  outputFormatGroup.add_argument("--sourcefile", "--file", "-f",
    action="store_true",
    help="print the file where triggering events are stored [%(default)s]"
    )
  outputFormatGroup.add_argument("--human", "--hr", "-H", action="store_true",
    help="print the event ID with extended labels [%(default)s]")
  outputFormatGroup.add_argument("--short", "-S", action="store_true",
    help="print the event ID with short labels [%(default)s]")

  generalOptGroup = parser.add_argument_group("General options")
  generalOptGroup.add_argument("--list", "-l", dest="doList", action="store_true",
    help="prints the available trigger data product tags, and exits")
  generalOptGroup.add_argument("--debug", action="store_true",
    help="enable additional diagnostic output")
  generalOptGroup.add_argument("--version", "-V", action="version",
    version=f"%(prog)s v{__version__} ({time.asctime(__date__)})",
    help="prints the version number"
    )
  
  args = parser.parse_args()
  
  logging.basicConfig(level=(logging.DEBUG if args.debug else logging.INFO))
  
  DefaultTriggerModuleLabel = args.deflabel
  # --- END ---- argument parsing  ---------------------------------------------
  
  sampleEvents = galleryUtils.makeEvent(args.InputFile)
  
  if args.doList: sys.exit(listTriggerDataProducts(sampleEvents))
  
  # process trigger specifications
  triggerSpecs = []
  if args.all or args.allToFiles:
    for tag in getAllTriggerTags(sampleEvents):
      triggerSpecs.append({
        'spec': tag,
        'destList': (buildDefaultOutputFile(tag) if args.allToFiles else None),
        })
    # for
  else:
    for spec in args.TriggerSpec:
      spec, sep, destList = spec.partition('->')
      if sep and not destList: destList = buildDefaultOutputFile(spec)
      triggerSpecs.append({ 'spec': spec, 'destList': destList, })
    # for
  # if all ... else
  
  sys.exit(findTriggeredEvents(sampleEvents, triggerSpecs, args))
  
# main
