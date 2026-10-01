# Generated SDK route is absent from shipped app

Product commit: 1ccbd0305cd9b639b06324fa136654d686483b7a.
Source: https://github.com/uditakankananonononono/sugarcode-ai

generate_sdk_client('crispr_opt', 'design_guides') emits POST /v1/modules/crispr_opt/design_guides. omega.api.app registers metadata, health and search GET routes, not a versioned module-operation POST route. In-process TestClient request to the generated route returns 404 twice. Positive control GET /modules/crispr_opt returns 200 and the correct slug. No external network POST was performed.

The existing Ecosystem suite checks generated code contains requests.post, an environment key and a string path. It does not execute the client against the shipped server. Thus its passing assertion is not an end-to-end SDK integration proof. The default domain is api.sugarcode.local, not a validated deployment.

Narrow finding: client/server contract mismatch in the shipped code. This report does not assert the state of an unseen deployment or recommend sending credentials to any URL. Product owner may implement an approval-safe operation route, generate in-process package calls instead, or explicitly label the client a future API template. No product changes here.
