<template>
	<div>
		<b-nav>
			<b-nav-item disabled><b>Quiz app</b></b-nav-item>
			<b-nav-item disabled>Counter: {{ numCorrect }}/ {{ numTotal }} </b-nav-item>
			<b-nav-item disabled>Field: {{ field }} </b-nav-item>
		</b-nav>
	</div>
</template>

<script>
export default {
	props: [
		'numCorrect',
		'numTotal',
		'field'
	]
}
</script>